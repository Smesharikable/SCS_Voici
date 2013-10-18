function [gmtest, input] = prepareTest()
    gmm = generateGMM(4, 2, 'diag', false, 10, 2);
    input = random(gmm, 100);
    figure;
    hold on;
    scatter(input(:, 1), input(:, 2)); 
    ezcontour(@(x, y)pdf(gmm, [x y]), [-15 15], [-15 15], 300)
    
    gmm = generateGMM(4, 2, 'diag', false, 10, 2);
    gmtest = gmmaptest(gmm);
end