function [] = plot_compared_errors(xVals, popModel, coeffs, addNoise, mulNoise, mulError)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%% Create the sample data set
numPoints = size(xVals, 2);

for i = 1:numPoints
    xVal = xVals(i);
    sampleData(i) = (poissrnd(popModel(coeffs, xVal)) + poissrnd(addNoise(xVal)))*normrnd(mulNoise(xVal), mulError(xVal));
end

% 'True Error'
% Variance at each point is just the point itself because this
% random variable follows a Poisson distribution
% Integrate the variance and take square root to get standard deviation
% of the integral

baseIntegral = trapz(xVals, popModel(coeffs, xVals));
addNoiseIntegral = trapz(xVals, addNoise(xVals));
mulNoiseIntegral = trapz(xVals, mulError(xVals).^2);
trueError = ((baseIntegral + addNoiseIntegral)*mulNoiseIntegral)^0.5;

maxy = 32;

figure;
subplot(1, 3, 1);
hold on;
plot(xVals, (popModel(coeffs, xVals)+addNoise(xVals)).*mulNoise(xVals), 'k', 'LineWidth', 4);
plot(xVals, sampleData, 'bo', 'LineWidth', 1.5);
ylim([0 maxy]);
hold off;
title(strcat('True Error:', {' '}, int2str(trueError)));
legend(gca, 'Population', 'Sample');

%% Fit using curve fit (fit_pulse)

[fitData, upperBound, lowerBound, cfError] = cf(xVals, sampleData, popModel, coeffs, 0);

subplot(1, 3, 2);
hold on;
plot(xVals, (popModel(coeffs, xVals)+addNoise(xVals)).*mulNoise(xVals), 'k', 'LineWidth', 4);
plot(xVals, fitData, 'm', 'LineWidth', 3);
plot(xVals, upperBound, 'r--', 'LineWidth', 3);
plot(xVals, lowerBound, 'b--', 'LineWidth', 3);
ylim([0 maxy]);
hold off;
title(strcat('Curve Fit Error:', {' '}, int2str(cfError)));
legend(gca, 'Population', 'Nonlinear Fit', 'Upper Bound', 'Lower Bound');

%% Fit using monte carlo

%% Need to include the error redistribution in the Monte-Carlo thing

[integrals, mcError] = mc(xVals, sampleData, popModel, coeffs);

subplot(1, 3, 3);
hold on;
histogram(integrals, 40);
hold off;
title(strcat('Monte-Carlo Error:', {' '}, int2str(mcError)));

end