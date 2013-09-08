% Calculation the probability of emission for a particular GMMS-model
function emis_pr = B_emis(numComp, GMMS, DATA)
fstMat = DATA{1};
countFtrComp = length(fstMat(:,1));
countGauss = GMMS.NComponents;
vecEmis = zeros(countFtrComp,1);
currEmis = zeros(countGauss,1);
emis_pr = 0;
matrix_feature = DATA{numComp};

for i = 1:countFtrComp
    for j = 1:countGauss
        currEmis(j,1) = GMMS.PComponents(j) * ...
        mvnpdf(matrix_feature(i,:), GMMS.mu(j,:), GMMS.Sigma(:,:,j));
    end
    for k = 1:countGauss
        vecEmis(i,1) = vecEmis(i,1) + currEmis(k,1);
    end
end
%  fprintf('%i', countGauss);
for i = 1:length(vecEmis)
    emis_pr = emis_pr + log(vecEmis(i,1));    
end    
end






