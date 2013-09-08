[fitgmm, input] = TryGMM(3);
%fitgmm.mu
gmm = fitwithmap(input, fitgmm, 16);
figure;
hold on;
scatter(input(:, 1), input(:, 2)) 
ezcontour(@(x, y)pdf(gmm, [x y]), [-15 15], [-15 15], 200)