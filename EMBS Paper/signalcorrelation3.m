clc
clear
close all

compmodel = readmatrix("SCG Displacement.csv");

t = compmodel(:,1);
disp = compmodel(:,2);

[dispfreq,dispamp] = FFT_fn(disp,74.074);

dispfilt = filter(filter_bandpass_disp(0.5,2.5),disp);
dispfilt = filter(filter_bandstop_2,dispfilt);

[peaks, x]  = findpeaks(dispfilt,74.074);


for i = 2:length(x)
    p(1) = 0;
    p(i) = x(i) - x(i-1);
end
meanp = mean(p);

bps = 1 / meanp;

bpm = bps * 60;


figure
plot(t,disp,"b")
hold on
scatter(x,peaks)
hold on
plot(t,dispfilt,"g")
title("ANSYS Displacement")
xlabel("Time (s)")
ylabel("Amplitude")
