function res = THD(Amplitudes,NumberHarmo)
    den = 0;
    num = 0;
    for i=1:NumberHarmo
        den = den + Amplitudes(i)^2;
    end
    num = den - Amplitudes(1)^2;
    
    res = sqrt(num)/sqrt(den)*100;
end