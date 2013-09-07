% UBM model creation
% Input dataMatrix - big amount of feature vectors
% for general model creation
function [UBM] = getUBM(dataMatrix, numComponents)
    options = statset('Display','final');
    UBM = gmdistribution.fit(dataMatrix, numComponents, 'CovType', 'diagonal', 'SharedCov', false,'Regularize', ...
        0.01,'Options', options);
end