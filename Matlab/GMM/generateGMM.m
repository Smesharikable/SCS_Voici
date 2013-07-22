function gmm = generateGMM(NComponents, NDimensions, CovType, isShared,...
        muBound, sigmaBound)
    
    mu = -muBound + 2 * muBound * rand(NComponents, NDimensions);
    
    if isShared
        K = 1;
    else 
        K = NComponents;
    end
    
    if strcmp(CovType, 'full')
        sigma = zeros(NDimensions, NDimensions, K);
        for i = 1:K
            sigma(:, :, i) = cov( -sigmaBound + 2 * sigmaBound * randn(NDimensions));
        end
    else
        sigma = 2 * sigmaBound * rand(1, NDimensions, K);
    end
    
    p = rand(1, NComponents);
    
    gmm = gmdistribution(mu, sigma, p);
end