              
function [summ] = denominator 
data = [2 2; 3 3];
NComponents = 2;
NDimensions = 2;
oldmu =[1 2; 2 3];
oldSigma = [2 0; 0 0.5];
compDensity = zeros(1, NComponents);
summ = zeros(1, NComponents);
compDensityMult = 1 / ((2 * pi)^(NDimensions/2) * (det(oldSigma)^(0.5)));
for t = 1:2
    for i = 1:NComponents
        for j = 1:NComponents;
            compDensity(j)= compDensityMult * exp(-0.5 * (data(t,:) - oldmu(j,:)) * inv(oldSigma) * (data(t,:) - oldmu(j,:))');
        end;
        summ(i) = sum(0.5 * compDensity(:));
        p(i) = (0.5 * compDensity(i)) / summ(i);
        Pr(t,i) = p(i);
    end
end