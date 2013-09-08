function [fitgmm, input] = TryGMM(count)
    % create gmm for input
    N = 4;
    mu = -10 + 20 * rand(N, 2);
    sigma = [2 0; 0 .5];
    % N Gaussian approximated by count Gaussian
    p = rand(1, N);
    gmm = gmdistribution(mu, sigma, p);
    input = random(gmm, 1000);

    %figure;
    %ezsurf(@(x, y)pdf(gmm, [x y]), [-15 15], [-15 15], 200)

    % applying EM algorithm
    mu = -10 + 20 * rand(N, 2);
    p = rand(1, N);
    fitgmm = gmdistribution(mu, sigma, p);
    %figure;
    %ezsurf(@(x, y)pdf(fitgmm, [x y]), [-15 15], [-15 15], 200)
end