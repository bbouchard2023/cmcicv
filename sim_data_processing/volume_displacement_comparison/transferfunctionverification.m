%==========================================================================
%% Transfer Function Verification
%
% Used to determine validity of coefficients for transfer function between
% LV volume and chest displacement
%
% Brendan Bouchard
% 20260330
% Last Updated: 20260330
%==========================================================================

clc
clear
close all

volfile = readmatrix("Ventricular_Volume.csv");
displfile = readmatrix("Displacement Data.csv");
accelfile = readmatrix("SCG Acceleration.csv");

t_vol = volfile(:,1);
t_vol = interp(t_vol,2);
t_vol = [0.0135; t_vol];
t_vol = t_vol(1:end-2);
vol = volfile(:,2);
vol = interp(vol,2);
vol = vol(1:end-1);
t_displ = displfile(:,1);
displ = displfile(:,4);
accel = accelfile(:,2);

% m     mass
% tau   propagation delay = 0
% c     damping coefficient
% k     stiffness



% gamma * deltaLV (ß) = ((2 .* displ) .* ((m./t_displ.^2) + (c./t_displ) + ((1/2).*k)));

%% Lungs

k_lungs = 1870; % average stiffness [Pa] (https://journals.physiology.org/doi/full/10.1152/ajplung.00415.2017#:~:text=After%20averaging%20all%20measurements%20within,location%20within%20the%20vascular%20tree.)
m_lungs = 0.565; % average lung mass in males [kg] (Normal organ weights in men: Part II—the brain, lungs, liver, spleen, and kidneys. American Journal of Forensic Medicine and Pathology,. https://doi.org/10.1097/PAF.0b013e31823d29ad)
c_crit_lungs = 2 * sqrt(k_lungs * m_lungs); % critical damping coefficient
damping_ratio_lungs = 2; % damping ratio (c/c_c) (https://journals.physiology.org/doi/abs/10.1152/ajplegacy.1956.186.1.142#:~:text=The%20natural%20frequency%20(9.6%2C%20S.E.,mechanics%20and%20ballistocardiography%20are%20discussed.)
c_lungs = damping_ratio_lungs * c_crit_lungs; % average damping coefficient [N*s/m]


Beta_lungs = ((2 .* displ) .* ((m_lungs./t_displ.^2) + (c_lungs./t_displ) + ((1/2).*k_lungs))); % gain * deltaLV



%% Ribcage


k_ribcage = 1870; % average stiffness [Pa] (https://pmc.ncbi.nlm.nih.gov/articles/PMC11835592/#:~:text=proapoptotic%20factor%20Bim.-,Results,1A).)
m_ribcage = 10.5; % average ribcage mass in males [kg] (https://orthoinfo.aaos.org/en/staying-healthy/healthy-bones-at-every-age/)
c_crit_ribcage = 2 * sqrt(k_ribcage * m_ribcage); % critical damping coefficient 
damping_ratio_ribcage = 0.089; % damping ratio (c/c_c) (https://www.tandfonline.com/doi/full/10.1080/13588265.2024.2348376#:~:text=Abstract,development%20of%20crash%20test%20dummies.)
c_ribcage = damping_ratio_ribcage * c_crit_ribcage; % average damping coefficient [N*s/m] 

Beta_ribcage = ((2 .* displ) .* ((m_ribcage./t_displ.^2) + (c_ribcage./t_displ) + ((1/2).*k_ribcage))); % gain * deltaLV


% %% Subcutaneous Fat
% 
% k_fat = 1870; % average stiffness [Pa] (https://pmc.ncbi.nlm.nih.gov/articles/PMC11835592/#:~:text=proapoptotic%20factor%20Bim.-,Results,1A).)
% m_fat = 10.5; % average ribcage mass in males [kg] (https://orthoinfo.aaos.org/en/staying-healthy/healthy-bones-at-every-age/)
% c_crit_fat = 2 * sqrt(k_fat * m_fat); % critical damping coefficient 
% damping_ratio_fat = 0.089; % damping ratio (c/c_c) (https://www.tandfonline.com/doi/full/10.1080/13588265.2024.2348376#:~:text=Abstract,development%20of%20crash%20test%20dummies.)
% c_fat = damping_ratio_fat * c_crit_fat; % average damping coefficient [N*s/m] 
% 
% Beta_fat = ((2 .* displ) .* ((m_fat./t_displ.^2) + (c_fat./t_displ) + ((1/2).*k_fat))); % gain * deltaLV
% 

Beta = Beta_lungs + Beta_ribcage;


f1 = figure;
subplot(2,1,1)
plot(t_displ,Beta_lungs,"m")
title("Transfer (Lungs)")
xlabel("Time (s)")
ylabel("Transfer (Γ*DeltaLV) [ß]")

subplot(2,1,2)
plot(t_displ,Beta_ribcage,"r")
title("Transfer (Ribcage)")
xlabel("Time (s)")
ylabel("Transfer (Γ*DeltaLV) [ß]")

saveas(f1,"IndividualTransfers.png")


f2 = figure;
plot(t_displ, Beta,"g")
title("Total Transfer")
xlabel("Time (s)")
ylabel("Transfer (Γ*DeltaLV) [ß]")

saveas(f2,"TotalTransfer.png")