function [freq,amp]=fft_fn(xn,fs)


Xk = fft(xn);
N = length(xn);
T = N/fs;

P2=abs(Xk); %two-seded spectrum amplitudes
P1=P2(1:N/2+1); % single-sided spectrum amplitudes
P1(2:end-1)=2*P1(2:end-1); % double the valuses between 0 and Fs/2
P1=P1/N; % divide by N
freq1=0:1/T:fs/2;  

figure;
stem(freq1,P1)
title('single sided spectrum')
xlabel(' Frequency Hz')
ylabel('2*|Xk|/N')

freq = freq1;   amp = P1;

end