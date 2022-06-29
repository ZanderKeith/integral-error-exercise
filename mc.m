function [integrals, mcError] = mc(xVals, data, popModel, coeffs)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

numRuns = 1000;

integrals = zeros(1, numRuns);

numPoints = size(xVals, 2);

standardDeviation = @(x) abs(x)^0.5;

for i = 1:numRuns
    newData = zeros(1, numPoints);
    for j = 1:numPoints
        %newData(j) = poissrnd(popModel(coeffs, data(j)));
        %newData(j) = poissrnd(data(j));
        %newData(j) = normrnd(data(j), standardDeviation(data(j)));
        
        %newData(j) = normrnd(data(j), standardDeviationScaled(data(j)));

        newData(j) = normrnd(data(j), standardDeviation(data(j)));
    end
    
    [fitData, ~, ~, ~] = cf(xVals, newData, popModel, coeffs, 1);

    integrals(i) = abs(trapz(xVals, fitData));
end

mcError = std(integrals);

end