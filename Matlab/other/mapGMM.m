%MAP
function mapGMM(NComponents, matrixType, SharedCov, mapIter, input1, input2)
    %Parametrs
    nObservations = 300;
    nCells = 200;
    %Random feature vectors
    input2 = randn(nObservations,2) * 7;
    options = statset('Display', 'final');
    %shows last iter statistic
    obj = gmdistribution.fit(input1, NComponents,'Start', 'randSample','CovType',matrixType,'SharedCov', SharedCov, 'Regularize',0.01, 'Options', options);
    gmmModel = gmdistribution(obj.mu, obj.Sigma, obj.PComponents);
    %InitFromInputData
    kComponents = gmmModel.NComponents;
    dDimension = gmmModel.NDimensions;
    muOld = gmmModel.mu;
    weightOld = gmmModel.PComponents;
    
    % create an array of covariation matrices
    if SharedCov == true
        sigmaOld = zeros(1, dDimension, kComponents);
        for i = 1:kComponents
            sigmaOld(:, :, i) = obj.Sigma;
        end
    else
        sigmaOld = obj.Sigma;
    end
    
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
    hold on;
    scatter(input1(:,1), input1(:,2));
 
    for k=1:mapIter
        for i=1:kComponents
            Pr = posterior(gmmModel, input2);
            for t=1:nObservations
                %Pr = weightOld(i)* mvnpdf(input2(t, :), muOld(i, :), sigmaOld(:,:,i)) / pdf(gmmModel, input2(t, :));
                n(i) = n(i) + Pr(t, i); 
                Ex(i, :) = Ex(i, :) + (input2(t, :) * Pr(t, i));
                Ex2(:, :, i) = Ex2(:, :, i) + (input2(t, :) .* input2(t,:))*Pr(t, i);
            end
            Ex(i, :) = Ex(i, :) / n(i);
            Ex2(:, :, i) = Ex2(:, :, i) / n(i);

            coef(i) = n(i) / (n(i) + 16);
            weights(i) = (coef(i) * n(i) / nObservations + (1 - coef(i)) * weightOld(i));
            means(i, :) = coef(i) * Ex(i) + (1 - coef(i)) * muOld(i, :);
            sigmas(:, :, i) = coef(i) * Ex2(:, :, i) + (1 - coef(i)) * (sigmaOld(:, :, i)...
                + (muOld(i, :) .* muOld(i, :))) - (Ex(i, :) .* Ex(i, :));
        end

        %Show MAPed gmmModel 
        gmmModel = gmdistribution(means, sigmaOld, weights);
        figure;
        ezcontour(@(x, y)pdf(gmmModel, [x y]), [-15 15], [-15 15], nCells);
        hold on;
        scatter(input2(:,1), input2(:,2));
        % next MAP iteration step
        muOld = means;
        %sigmaOld = sigmas;
        weightOld = weights;
        n = zeros(1,kComponents);
        Ex = zeros(kComponents, dDimension);
        Ex2= zeros(1, dDimension, kComponents);
    end    
end


    
