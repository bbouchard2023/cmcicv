%==========================================================================
%% Signal Comparison
%
% Compares displacement between computational model and SCG data, and
% evaluates correlation between ventricular volume and chest displacement
% 
%
% Brendan Bouchard
% 20260122
% Last Updated: 20260122
%==========================================================================

clc
clear
close all

%% Model/SCG Displacement Comparison

ansys_displacement = readmatrix("displacement_data.xlsx");
ansys_t = ansys_displacement(:,1);
ansys_disp = ansys_displacement(:,4);

fig = figure;
fig.WindowState = "maximized";
plot(ansys_t,ansys_disp,"m")
title("Computational Model/SCG Displacement Comparison")
xlabel("Time (s)")
ylabel("Displacement (m)")
legend("Computational Model Displacement")


%% Ventricular Volume/Chest Displacement Correlation

files = dir("volume*.mat");
frameNum = arrayfun(@(f) sscanf(f.name,'volume%d.mat'), files);
[frameNum, idx] = sort(frameNum);
files = files(idx);
numfiles = numel(files);
vol = zeros(numfiles,1);
dt = 0.034;

for l = 1:numfiles
    S = load(files(l).name);
    J = S.J;
    vol(l) = nnz(J);
end

t = (frameNum - frameNum(1)) * dt;

fig2 = figure;
fig2.WindowState = "maximized";
subplot(2,1,1)
plot(t,vol,'m');
title("Cardiac Volume")
xlabel("Time (s)")
ylabel("Volume (mL)")
xlim([t(1) ansys_t(end)])
subplot(2,1,2)
plot(ansys_t,ansys_disp,"g")
title("Computational Model Displacement")
xlabel("Time (s)")
ylabel("Displacement (m)")

