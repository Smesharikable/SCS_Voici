% input - features vector
% gausscount - pozitive integer of gmm components
% srcgmm - source gaussian mixture model
% relfactor - pozitive integer representing measure of influence of new
%   data to adaptation source model 
function gmm = fitwithmap(srcgmm, relfactor)
    input = random(srcgmm, 1000);
    if strcmp(srcgmm.CovType, 'full')
        gmm = gmmMap(input, srcgmm, relfactor, @diag);
    else 
        gmm = gmmMap(input, srcgmm, relfactor, @emptyFunc);
        input
    end
    
end

function gmm = gmmMap(input, srcgmm, relfactor, correctDimension)
    mu = srcgmm.mu;
    w = srcgmm.PComponents;
    N = srcgmm.NComponents;
    D = srcgmm.NDimensions;
    T = length(input);
    CovDim = length(srcgmm.Sigma(:, 1, 1));
    
    % create an array of covariation matrices
    if srcgmm.SharedCov == true
        sigma = zeros(CovDim, D, N);
        for i = 1:N
            sigma(:, :, i) = srcgmm.Sigma;
        end
    else
        sigma = srcgmm.Sigma;
    end
    
    weights = zeros(1, N);
    sigmas = zeros(CovDim, D, N);
    means = zeros(N, D);
    coef = zeros(1, N);
    
    % compute nessesary statistics
    for i = 1:N
        for t = 1:T
            temp = w(i) * mvnpdf(input(t, :), mu(i, :), sigma(:, :, i)) / pdf(srcgmm, input(t, :));
            weights(i) = weights(i) + temp;
            temp = input(t, :) .* temp;
            means(i, :) = means(i, :) + temp;
            sigmas(:, :, i) = sigmas(:, :, i) + correctDimension(input(t, :) .* temp);
        end
        means(i, :) = means(i, :) ./ weights(i);
        sigmas(:, :, i) = sigmas(:, :, i) ./ weights(i);
        % parameters adaptation
        coef(i) = weights(i) / (weights(i) + relfactor);
        weights(i) = (coef(i) * weights(i) / T + (1 - coef(i)) * w(i));
        means(i, :) = coef(i) .* means(i) + (1 - coef(i)) .* mu(i, :);
        sigmas(:, :, i) = coef(i) .* sigmas(:, :, i) + (1 - coef(i)) .* (sigma(:, :, i)...
            + correctDimension(mu(i, :) .* mu(i, :))) - correctD  imension(means(i, :) .* means(i, :));
    end
           
    gmm = gmdistribution(means, sigmas, weights);
    figure;
    ezsurf(@(x, y)pdf(gmm, [x y]), [-15 15], [-15 15], 200);
    
end

function output = emptyFunc(input)
    output = input;
end