%==========================================================================
%% Cardiovascular Volume/Displacement Comparison
%
% Compares cardiovascular volume and displacement
%
% Brendan Bouchard
% 20260131
% Last Updated: 20260204
%==========================================================================

clc
clear
close all




dispdata = readmatrix("Displacement Data.csv");
displ = dispdata(:,2);
t = dispdata(:,1);


voldata = readmatrix("Ventricular Volume.xlsx");
vol = voldata(:,2);
vol = [vol; vol; vol; vol; vol];
vol = interp(vol,2);
vol = vol(2:end);

rho = corrcoef(vol,displ);
rho = rho(1,2);

disp(rho)




[volzeros, tvol] = findpeaks(vol);
[displzeros, tdispl] = findpeaks(displ);


test = (volzeros > 154400) | (volzeros < 150000 & volzeros > 140000);
volzeros = volzeros(test);
tvol = tvol(test);

voldisplag = abs(t(tvol) - t(tdispl));
avglag = mean(voldisplag);

disp(avglag)

fig = figure;
fig.Name = "Ventricular Volume and Chest Displacement Comparison";

title("Ventricular Volume and Chest Displacement Comparison")
yyaxis left
plot(t,displ,"m")

xlabel("Time (s)")
ylabel("Displacement (mm)")
hold on
scatter(t(tdispl),displzeros)
yyaxis right
plot(t,vol,"g")
ylabel("Volume (mm^3)")
scatter(t(tvol),volzeros,"yx")


legend("Displacement","Cardiac Volume","Volume Peaks","Displacement Peaks")
text(4,50000,"\rho = " + num2str(rho))


