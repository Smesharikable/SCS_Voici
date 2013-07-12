%MAP
function mapGMM(NComponents)
    %Parametrs
    nObservations = 1000;
    nCells = 200;
    %Random feature vectors
    input = rand(nObservations,2) * 20;
    options = statset('Display', 'final');
    %shows last iter statistic
    obj = gmdistribution.fit(input, NComponents,'Start', 'randSample','CovType','diagonal','SharedCov', false, 'Regularize',0.01, 'Options', options);
    gmmModel = gmdistribution(obj.mu, obj.Sigma, obj.PComponents);
    %InitFromInputData
    kComponents = gmmModel.NComponents;
    dDimension = gmmModel.NDimensions;
    muOld = gmmModel.mu;
    sigmaOld = gmmModel.Sigma;
    weightOld = gmmModel.PComponents;
    
    %Declaire empty variables
    n = zeros(1,kComponents);
    Ex = zeros(kComponents, dDimension);
    Ex2= zeros(1, dDimension, kComponents);
    
    %Declaire empty parametrs for MAP GMM
    weights = zeros(1, kComponents);
    sigmas = zeros(1, dDimension, kComponents);
    means = zeros(kComponents, dDimension);
    coef = zeros(1, kComponents);
    
    %Show start gmmModel
    figure;
    ezcontour(@(x, y)pdf(gmmModel, [x y]), [-15 15], [-15 15], nCells);
    
    %Find coef
    for i=1:kComponents
        for t=1:nObservations
            Pr = weightOld(i)* mvnpdf(input(t, :), muOld(i, :), sigmaOld(:,:,i)) / pdf(gmmModel, input(t, :));
            n(i) = n(i) + Pr; 
            Ex(i, :) = Ex(i, :) + (input(t, :) * Pr);
            Ex2(:, :, i) = Ex2(:, :, i) + (input(t, :) .* input(t,:))*Pr;
        end
        Ex(i, :) = Ex(i, :) / n(i);
        Ex2(:, :, i) = Ex2(:, :, i) / n(i);
        
        coef(i) = n(i) / (n(i) + 16);
        weights(i) = (coef(i) * n(i) / nObservations + (1 - coef(i)) * weightOld(i));
        means(i, :) = coef(i) * Ex(i) + (1 - coef(i)) * muOld(i, :);
        sigmas(:, :, i) = coef(i) .* Ex2(:, :, i) + (1 - coef(i)) .* (sigmaOld(:, :, i)...
            + (muOld(i, :) .* muOld(i, :))) - (Ex(i, :) .* Ex(i, :));
    end
    
    %Show MAPed gmmModel with help coef from imaginary UBM
    GMMLastIter = gmdistribution(means, sigmas, weights);
    figure;
    ezcontour(@(x, y)pdf(GMMLastIter, [x y]), [-15 15], [-15 15], nCells);
    
    %Change values | update data | more MAP iter
    

end
    
