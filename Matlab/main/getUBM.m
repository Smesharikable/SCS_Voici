% UBM model creation
% Input dataMatrix - big amount of feature vectors
% which helps to create an general model
function [ubm] = getUBM(dataMatrix, numComponents)
    options = statset('Display','final');
    ubm = gmdistribution.fit(dataMatrix, numComponents, 'CovType', 'diagonal', 'SharedCov', false,'Regularize', ...
        0.01,'Options', options);
end