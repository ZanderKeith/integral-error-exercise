close all
rng('default') % for reproducibility

xVals = 0:200;

coeffs = 40;
constantPopulationFunction = @(g, x) g(1) + (x*0);

p = Population(xVals, constantPopulationFunction, coeffs);
p.triplePlot('Constant Signal');

