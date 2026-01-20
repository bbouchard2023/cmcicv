clc
clear
close all

compmodel = readmatrix("SCG Displacement.csv");

t = compmodel(:,1);
disp = compmodel(:,2);

[dispfreq,dispamp] = FFT_fn(disp,74.074);

dispfilt = filter(filter_bandpass_disp,disp);



figure
plot(t,disp,"b")
hold on
plot(t,dispfilt,"g")
title("ANSYS Displacement")
xlabel("Time (s)")
ylabel("Amplitude")
legend("Raw Displacement", "Bandpass Filtered Displacement")
hold off
