function [input1, input2] = gmfit
    % Parametrs
    kComponents = 2;
    nObservations = 300;
    dDimension = 2;
    nCells = 300;
    % Init Step
    % for [-10;10] diapason
    mu = -10 + rand(kComponents,dDimension) * 20;
    % different diagonal Cov matrix per Gaussian
    sigma = cat(3,[2 0;0 .5],[1 0;0 1]);
    p = rand(1,kComponents);
    GMM = gmdistribution(mu, sigma, p);
    input1 = random(GMM, nObservations); 
    
    %figure;
    %ezsurf(@(x, y)pdf(GMM, [x y]), [-15 15], [-15 15], nCells);
    
    % EM - algorithm
    options = statset('Display', 'final'); %shows last iter statistic
    obj = gmdistribution.fit(input1,kComponents,'Start', 'randSample','CovType','diagonal','SharedCov', false, 'Regularize',0.01, 'Options', options);
    ResultGMM = gmdistribution(obj.mu, obj.Sigma, obj.PComponents);
    input2 = random(ResultGMM, nObservations);
end





