function [fitgmm, input] = TryGMM(count)
    % create gmm for input
    N = 4;
    mu = -10 + 20 * rand(N, 2);
    sigma = [2 0; 0 .5];
    % N Gaussian approximated by count Gaussian
    p = rand(1, N);
    gmm = gmdistribution(mu, sigma, p);

    figure;
    ezcontour(@(x, y)pdf(gmm, [x y]), [-15 15], [-15 15], 200)

    % applying EM algorithm
    input = random(gmm, 1000);
    fitgmm = gmdistribution.fit(input, count, 'Start', 'randSample', 'CovType', 'diagonal','SharedCov', false, 'Regularize', 0.01);
    figure;
    ezcontour(@(x, y)pdf(fitgmm, [x y]), [-15 15], [-15 15], 200)
end