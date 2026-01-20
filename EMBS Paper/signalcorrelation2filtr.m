clc
clear
close all

compmodel = readmatrix("SCG Displacement.csv");

t = compmodel(:,1);
disp = compmodel(:,2);
disp1 = disp;
[dispfreq,dispamp] = FFT_fn(disp,74.074);

dispfilt11 = filter(filter_lowpass_20,disp1);

for low = 0.4:0.1:1.3
    for high = 3.5:0.1:4.4
        dispfilt = filter(filter_bandpass_disp(low,high),disp);
        delta = abs((disp - dispfilt)./disp);
    end
 
end

maxcorr = min(delta);


figure;
plot(t,disp,"b")
hold on
plot(t,dispfilt11,"g")
plot(t,dispfilt,"r")
l = "ANSYS Displacement, Low: " + num2str(low) + " , High: " + num2str(high) ;
title(l)
xlabel("Time (s)")
ylabel("Amplitude")
legend("Raw Displacement", "lowpass","bandpass")
hold off

