function Mmax = myMax(M)
Mmax = -1000000;
for i = 1:length(M)
    if M(i,1) > Mmax && M(i,1) ~= 0
        Mmax = M(i,1);
    end
end
if Mmax == -1000000
        Mmax = 0;
end
end