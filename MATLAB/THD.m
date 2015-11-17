function res = THD(A,N)
    den = 0;
    for i=1:N
        den = den + A(i)^2;
    end
    num = den - A(1)^2;
    
    res = sqrt(num)/sqrt(den)*100;
end