function [ res ] = amplituder(x, TIME_SIM, fr)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
NFFT = length(x);
%f = Fs*(0:NFFT/2)/NFFT;
X = fft(x,NFFT)/(NFFT);
X=2*abs(X(1:NFFT/2+1)); 
res=X(fr*(TIME_SIM-1)+1);
end

