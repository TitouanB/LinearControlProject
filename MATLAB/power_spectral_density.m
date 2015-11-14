function [Pxx,f] = power_spectral_density(x,Fs)
    % Generate the frequency vector
    NFFT = length(x);
    f = Fs*(0:NFFT/2)/NFFT;

    % Calculate the one sided power spectral density
    X = fft(x,NFFT);
    Px = abs(X).^2/(Fs*NFFT);
    Pxx = Px(1:NFFT/2+1);
    Pxx(2:end-1) = 2*Pxx(2:end-1);

end