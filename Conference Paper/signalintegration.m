clc
clear

fid = readtable("DirectionalAccel.xlsx");
arr = table2array(fid);

acc = arr(:,3);
accfilt = filter(filter_lowpass_08,acc);
t = arr(:,1);
fs = 74.074;

dzdt = cumtrapz(accfilt,t);

z = cumtrapz(dzdt,t);

zfilter = filter(filter_lowpass_08,z);

plot(t,zfilter,"m");
hold on
plot(t,accfilt,"g");
xlabel("Time (s)")
ylabel("Amplitude")
legend("Twice Integrated Position","Acceleration")
title("Twice Integrated Position from Acceleration (SCG Computational Model)")
writematrix(z,"IntegratedPosition.csv")
