[fitgmm, input] = TryGMM(3);
%fitgmm.mu
gmm = fitwithmap(input, fitgmm, 16);
figure;
ezsurf(@(x, y)pdf(gmm, [x y]), [-15 15], [-15 15], 200)