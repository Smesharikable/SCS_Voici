function Imax = argmax(M)
Mmax = -1000000;
Imax = 0;
for i = 1:length(M)
    if M(i,1) > Mmax && M(i,1) ~= 0
        Mmax = M(i,1);
        Imax = i;
    end
end
end