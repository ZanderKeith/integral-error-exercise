close all
rng('default') % for reproducibility

xVals = 0:200;

%% Noise sources
photonErrorFunc = @(g, x) g(1) + (x*0);
photonErrorCoeffs = 15;

gainFunc = @(x) 1 + (x*0);
gainSigma = 0.10;

electronErrorFunc = @(g, x) g(1) + (x*0);
electronErrorCoeffs = 20;

%% Constant population function
coeffs = 40;
constantPopulationFunction = @(g, x) g(1) + (x*0);

% % No noise
% 
% c0 = Population(xVals, constantPopulationFunction, coeffs);
% c0.triplePlot('Constant Signal');
% 
% % add photon noise
% 
% c1 = Population(xVals, constantPopulationFunction, coeffs, ...
%     photonErrorFunc, photonErrorCoeffs);
% c1.triplePlot('Constant Signal');

% c2 = Population(xVals, constantPopulationFunction, coeffs, ...
%     photonErrorFunc, photonErrorCoeffs, ...
%     gainFunc, gainSigma);
% c2.triplePlot('Constant Signal');

c3 = Population(xVals, constantPopulationFunction, coeffs, ...
    photonErrorFunc, photonErrorCoeffs, ...
    gainFunc, gainSigma, electronErrorFunc, electronErrorCoeffs);
c3.triplePlot('Constant Signal');

%% Linear population function

coeffs = [0.4 20];
linearPopulationFunction = @(g, x) g(1)*x + g(2);

l3 = Population(xVals, linearPopulationFunction, coeffs, ...
    photonErrorFunc, photonErrorCoeffs, ...
    gainFunc, gainSigma, electronErrorFunc, electronErrorCoeffs);
l3.triplePlot('Linear Signal');

%% Gaussian population function

coeffs = [40 100 50 0];
gaussianPopulationFunction = @(g, x) g(1).*exp(-((x - g(2))/g(3)).^2) + g(4);

g3 = Population(xVals, gaussianPopulationFunction, coeffs, ...
    photonErrorFunc, photonErrorCoeffs, ...
    gainFunc, gainSigma, electronErrorFunc, electronErrorCoeffs);
g3.triplePlot('Gaussian Signal');