function [fitgmm, input] = TryGMM(count)
    % create gmm for input
    mu = [-1, -10; 3, -5; 4, 9; -3, 5];
    sigma = [2 0; 0 .5];
    p = rand(1, 4);
    gmm = gmdistribution(mu, sigma, p);

    figure;
    ezsurf(@(x, y)pdf(gmm, [x y]), [-15 15], [-15 15], 200)

    % applying EM algorithm
    input = random(gmm, 1000);
    fitgmm = gmdistribution.fit(input, count, 'Start', 'randSample', 'CovType', 'diagonal','SharedCov', true, 'Regularize', 0.01);
    figure;
    ezsurf(@(x, y)pdf(fitgmm, [x y]), [-15 15], [-15 15], 200)
end