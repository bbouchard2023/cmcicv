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
delta = abs(vol - displ);

rho = corrcoef(vol,displ);
rho = rho(1,2);

disp(rho)

fig = figure;
fig.Name = "Ventricular Volume and Chest Displacement Comparison";

title("Ventricular Volume and Chest Displacement Comparison")
yyaxis left
plot(t,displ,"m")
xlabel("Time (s)")
ylabel("Displacement (mm)")
hold on
yyaxis right
plot(t,vol,"g")
ylabel("Volume (mm^3)")
legend("Displacement","Cardiac Volume")
text(4,50000,"\rho = " + num2str(rho))

