classdef Population < handle
    %Population Represents a random variable as well as noise sources which may affect it
    %   Detailed explanation goes here

    properties
        % X values
        xVals;

        % Population function 
        % Assumed to be a poisson random variable
        popFunction;
        popCoeffs;
        
        % Phton error (added directly to population)
        photonErrorFunc;
        photonErrorCoeffs;
        
        % Gain function (Multiply population * photon error)
        gainFunc;
        gainSigma;
        
        % Electronics error (added after gain multiplication)
        electronErrorFunc;
        electronErrorCoeffs;

        % 'Population sample' used by all graphs
        popSample;
    end

    methods
        %% Constructor
        function obj = Population(xVals, popFunction, popCoeffs, ...
                photonErrorFunc, photonErrorCoeffs, ...
                gainFunc, gainSigma, ...
                electronErrorFunc, electronErrorCoeffs)
            %Population Construct an instance of this class
            %   Detailed explanation goes here
            obj.xVals = xVals;

            obj.popFunction = popFunction;
            obj.popCoeffs = popCoeffs;
            
            if nargin <= 3
                % No function input, set photon error to 0 everywhere
                obj.photonErrorFunc = @(g, x) g(1)*0 + (x*0);
                obj.photonErrorCoeffs = 0;
            else
               obj.photonErrorFunc = photonErrorFunc;
               obj.photonErrorCoeffs = photonErrorCoeffs;
            end
            
            if nargin <= 5
                % No function input, set gain to 1 everywhere
                obj.gainFunc = @(x) 1 + (x*0);
                obj.gainSigma = 0;
            else
                obj.gainFunc = gainFunc;
                obj.gainSigma = gainSigma;
            end
            
            if nargin <= 7
                obj.electronErrorFunc = @(g, x) g(1)*0 + (x*0);
                obj.electronErrorCoeffs = 0;
            else
                obj.electronErrorFunc = electronErrorFunc;
                obj.electronErrorCoeffs = electronErrorCoeffs;
            end
            
        end

        %% Helper functions
        function sample = makeSample(obj)
            % Make a sample of the population with all noise sources included
            population = poissrnd(obj.popFunction(obj.popCoeffs, obj.xVals));
            photonNoise = poissrnd(obj.photonErrorFunc(obj.photonErrorCoeffs, obj.xVals));
            gain = normrnd(obj.gainFunc(obj.xVals), obj.gainSigma);
            electronNoise = poissrnd(obj.electronErrorFunc(obj.electronErrorCoeffs, obj.xVals));

            sample = ((population + photonNoise).*gain) + electronNoise;
        end

        function [trueError] = findTrueError(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here

            popVariance = @(x) obj.popFunction(obj.popCoeffs, x);
            photonVariance = @(x) popVariance(x) + obj.photonErrorFunc(obj.photonErrorCoeffs, x);
            photonMean = @(x) photonVariance(x); % Poisson variables

            % This gets a little kooky but I think I got it. TODO comment what's going on here
            gainVariance = obj.gainSigma^2;
            gainMean = 1; % ASSUMING gain mean stays 1 because we only really care about the variation, so variation around unity is fine

            productVariance = @(x) ((photonVariance(x) + (photonMean(x).^2)).*(gainVariance+gainMean)) - (photonMean(x).^2);

            electronVariance = @(x) productVariance(x) + obj.electronErrorFunc(obj.electronErrorCoeffs, x);

            varianceIntegral = trapz(obj.xVals, electronVariance(obj.xVals));

            trueError = varianceIntegral^0.5;
        end

        function [expectedValue] = getExpectedValue(obj)
            x = obj.xVals;

            pop = obj.popFunction(obj.popCoeffs, x);
            photonNoise = obj.photonErrorFunc(obj.photonErrorCoeffs, x);
            gain = obj.gainFunc(x);
            electronNoise = obj.electronErrorFunc(obj.electronErrorCoeffs, x);

            expectedValue = ((pop+photonNoise).*gain) + electronNoise;
        end

        function [integrals, monteCarloError] = getMonteCarloError(obj)
            numRuns = 10000;
            
            integrals = zeros(1, numRuns);

            for i = 1:numRuns
                newData = makeSample(obj);
                [fitData, ~, ~, ~] = cf(obj.xVals, newData, obj.popFunction, obj.popCoeffs, 1);

                integrals(i) = abs(trapz(obj.xVals, fitData));
            end
            monteCarloError = std(integrals);

        end

        %% Plotting

        function [] = triplePlot(obj, superTitle)
            figure;
            
            obj.plotTrueError();
            obj.plotCurveFitError();
            obj.plotMonteCarloDistribution();

            photonString = ['Stray Light: ', int2str(obj.photonErrorCoeffs(1))];
            gainString = ['Gain Sigma: ', int2str((obj.gainSigma)*100), '%'];
            electroString = ['Electronics Noise: ', int2str(obj.electronErrorCoeffs(1))];

            sgtitle({superTitle, photonString, gainString, electroString});
        end

        function [] = plotTrueError(obj)
            % Create a sample data set
            obj.popSample = obj.makeSample();
            
            x = obj.xVals;
            
            y = obj.getExpectedValue();
            
            subplot(1, 3, 1);
            hold on;
            plot(x, y, 'k', 'LineWidth', 4);
            plot(x, obj.popSample, 'bo', 'LineWidth', 1.5);
            ylim([0 max(obj.popSample)]);
            hold off;
            title(strcat('True Error:', {' '}, int2str(obj.findTrueError())));
            legend(gca, 'Expected Value', 'Sample');
        end % end plotTrueError

        function [] = plotCurveFitError(obj)
            x = obj.xVals;

            [fitData, upperBound, lowerBound, cfError] = cf(x, obj.popSample, obj.popFunction, obj.popCoeffs, 0);
            
            y = obj.getExpectedValue();

            subplot(1, 3, 2);
            hold on;
            plot(x, y, 'k', 'LineWidth', 4);
            plot(x, fitData, 'm', 'LineWidth', 3);
            plot(x, upperBound, 'r--', 'LineWidth', 3);
            plot(x, lowerBound, 'b--', 'LineWidth', 3);
            ylim([0 max(obj.popSample)]);
            hold off;
            title(strcat('Curve Fit Error:', {' '}, int2str(cfError)));
            legend(gca, 'Expected Value', 'Nonlinear Fit', 'Upper Bound', 'Lower Bound');


        end % end plotCurveFitError

        function [] = plotMonteCarloDistribution(obj)
            x = obj.xVals;

            [integrals, mcError] = obj.getMonteCarloError();

            subplot(1, 3, 3);
            hold on;
            histogram(integrals, 20);
            hold off;
            title(strcat('Monte-Carlo Error:', {' '}, int2str(mcError)));
            
        end
    end
end