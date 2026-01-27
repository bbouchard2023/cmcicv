%==========================================================================
%% 3D Volume Interpolator
%
% Interpolate 3D volume data in Z direction
%
% Original: Peshala T Gamage
% Refactored: Brendan Bouchard
% 20251125
% Last Updated: 20251125
%==========================================================================

clc
clear



for fr = 2:30 % for frames 2 - 30
    
    Vname = ['volume' num2str(fr) '.mat']; % defines file name for volume
    Vol = load(Vname); % loads volume .mat file
    
    J = Vol.J; % calls a cell array as variable J
    
    zold = linspace(1,  size(J,2), size(J,2))'; % linspace of original size of data
    znew = linspace(1 , size(J,2), size(J,2) * 4)'; % linspace for interpolated data size (4x the original number of points)

    array = reshape(J, [], 14).'; % automatically calculates the length of the second dimension
    new = interp1(zold, array, znew).'; % interpolates between real values of J
    J_new = reshape(new, [256,256,56]); % puts interpolated J into a new array
    clear J % clears original J array from workspace
    
    Vname_new = ['volumeinterp' num2str(fr) '.mat']; % names new .mat file with interpolated J array
    save(Vname_new,'J_new'); % saves new .mat file
end