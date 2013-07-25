function [ubm] = getUBM(dataMatrix, numComponents)
    options = statset('Display','final');
    ubm = gmdistribution.fit(dataMatrix, numComponents, 'CovType', 'diagonal', 'SharedCov', false,'Regularize', ...
        0.0001,'Options', options);
end