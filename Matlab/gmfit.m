function [ResultGMM, input] = gmfit(kComponents, approxComponents)
    % Parametrs
    nObservations = 1000;
    dDimension = 2;
    nCells = 300;
    % Init Step
    % for [-10;10] diapason
    mu = -10 + rand(kComponents,dDimension) * 20;
    % different diagonal Cov matrix per Gaussian
    sigma = zeros(1, dDimension, kComponents);
    for i=1:kComponents
        sigma(:, :, i) = diag(cov(-1 + 2 * randn(dDimension)))';
    end
    p = rand(1,kComponents);
    GMM = gmdistribution(mu, sigma, p);
    input = random(GMM, nObservations); 
    
    %figure;
    %ezsurf(@(x, y)pdf(GMM, [x y]), [-15 15], [-15 15], nCells);
    
    % EM - algorithm
    options = statset('Display', 'final'); %shows last iter statistic
    obj = gmdistribution.fit(input,approxComponents,'Start', 'randSample','CovType','diagonal','SharedCov', false, 'Regularize',0.01, 'Options', options);
    ResultGMM = gmdistribution(obj.mu, obj.Sigma, obj.PComponents);
    
    %figure;
    %ezsurf(@(x, y)pdf(ResultGMM, [x y]), [-15 15], [-15 15], nCells);
end





