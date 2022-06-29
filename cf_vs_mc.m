%% Curve fit error vs Monte-Carlo error sandbox because this has been bugging me for far too long
close all
rng('default') % for reproducibility

%% Set up vars
numPoints = 100;
xVals = 1:numPoints;

addNoiseFunction = @(x) 8 + (x*0); % additive error
mulNoiseFunction = @(x) 1.5 + (x*0); % multiplicative error
mulNoiseErrorFunction = @(x) 0.2 + (x*0);

%% Constant model
coeffs = 8;
constantPopulationFunction = @(g, x) g(1) + (x*0);
plot_compared_errors(xVals, constantPopulationFunction, coeffs, addNoiseFunction, mulNoiseFunction, mulNoiseErrorFunction);

%% Linear model
coeffs = [0.05 6];
linearPopulationFunction = @(g, x) g(2) + (x * g(1));
plot_compared_errors(xVals, linearPopulationFunction, coeffs, addNoiseFunction, mulNoiseFunction, mulNoiseErrorFunction);

%% Gaussian model
coeffs = [12 50 20];
coeffs = [12 50 20 2];
cubicPopulationFunction = @(g, x) g(1).*exp(-((x - g(2))/g(3)).^2) + g(4);
plot_compared_errors(xVals, cubicPopulationFunction, coeffs, addNoiseFunction, mulNoiseFunction, mulNoiseErrorFunction);