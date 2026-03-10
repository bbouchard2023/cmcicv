%==========================================================================
%% Integration Comparison
%
% Compares Twice Integrated Acceleration with Displacement from Ansys
% Simulation
%
% Brendan Bouchard
% 20260224
% Last Updated: 20260224
%==========================================================================

clc
clear
close all


displacement_table = readmatrix("Point Displacement.csv"); % imports displacement data from file
acceleration_table = readmatrix("Point SCG.csv"); % imports acceleration data from file

t_displacement = displacement_table(:,1); % time vector for displacement
displacement = displacement_table(:,2); % displacement vector

t_acceleration = acceleration_table(:,1); % time vector for acceleration
acceleration = acceleration_table(:,2); % acceleration vector


velocity = cumtrapz(acceleration); % velocity vector derived by numerically integrating acceleration from Ansys simulation

int_displacement = cumtrapz(velocity); % displacement vector derived by numerically integrating acceleration from Ansys simulation twice

f = figure;
f.WindowState = "maximized";

subplot(4,1,1)
plot(t_acceleration, acceleration, "g")
title("Acceleration from Ansys Simulation")
xlabel("Time (s)")
ylabel("Acceleration (mm/s^2)")

subplot(4,1,2)
plot(t_acceleration, velocity, "c")
title("Velocity Integrated from Ansys Simulation Acceleration Vector")
xlabel("Time (s)")
ylabel("Velocity (mm/s)")

subplot(4,1,3)
plot(t_displacement,displacement,"r")
title("Displacement from Ansys Simulation")
xlabel("Time (s)")
ylabel("Displacement (mm)")

subplot(4,1,4)
yyaxis left
plot(t_displacement,displacement,"r")
ylabel("Displacement (mm)")

hold on

yyaxis right
plot(t_acceleration,int_displacement,"b")
xlabel("Time (s)")
ylabel("Displacement (mm)")
title("Comparison Between Twice Integrated Acceleration (Displacement) and Displacement from Ansys Simulation")
legend("Displacement from Ansys Simulation","Twice Integrated Acceleration (Displacement)")

saveas(f,"IntegrationComparison.png")

Fs = 74.074;
[displfreq, displamp] = fft_fn(displacement,74.074);


