function [fitData, upperBound, lowerBound, cfError] = cf(xVals, data, popModel, guess, skip)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[coefficientValues, R, J, CovB, MSE, ErrorModelInfo] = nlinfit(xVals, data, popModel, guess);
% Get fit and 68% (1 sigma) confidence intervals
[fitData, delta] = nlpredci(popModel, xVals, coefficientValues, R, 'Covar', CovB, 'MSE', MSE, 'SimOpt', 'on', 'Alpha', 1 - 0.68);

if skip == 0
    upperBound = fitData + delta;
    lowerBound = fitData - delta;
    
    % Integrate upper and lower bounds, take difference, divide by 2
    upperVal = trapz(xVals, upperBound);
    lowerVal = trapz(xVals, lowerBound);
    
    % Error produced by confidence interval method is proportional to degrees of freedom
    cfError = (upperVal - lowerVal) / (2 * size(guess, 2));

else 
    upperBound = 0;
    lowerBound = 0;
    cfError = 0;
end% end skip

end