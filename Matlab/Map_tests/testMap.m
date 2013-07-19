gmm = generateGMM(4, 2, 'diag', false, 10, 2);
input = random(gmm, 100);
figure;
hold on;
scatter(input(:, 1), input(:, 2)) 
ezcontour(@(x, y)pdf(gmm, [x y]), [-15 15], [-15 15], 200)

N = 20;
statistic = zeros(1, N);
for i = 1:N
    gmm = generateGMM(4, 2, 'diag', false, 10, 2);
    test = gmmaptest(gmm);
    %test.visualizing = true;

    t1 = test.fitbymeans(input, 3);
    t2 = test.fitbymeans2(input, 3);

    %t1 = test.fitbyweight(input, 3);
    %t2 = test.fitbyweight2(input, 3);

    y1 = pdf(t1.gmmout, input);
    y2 = pdf(t2.gmmout, input);
    y = sum(y1);
    statistic(i) = (y - sum(y2)) / y;
end

sum(statistic) / N
    