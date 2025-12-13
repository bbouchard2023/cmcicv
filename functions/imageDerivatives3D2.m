%==========================================================================
%% 3D Derivative Between 2 3D Images
%
% Computes 3D derivatives between two 3D images
%
% Original: Peshala T Gamage
% Refactored: Brendan Bouchard
% 20251125
% Last Updated: 20251125
%==========================================================================



function [Ix, Iy, Iz, It] = imageDerivatives3D2(image1, image2, image3, delt)


dx = zeros(2,2,2); % preallocates dx array
dx(:,:,1) = [-1 1; -1 1 ]; % 
dx(:,:,2) = [-1 1; -1 1 ]; % 
dx = 0.25 * dx; % 

dy = zeros(2,2,2); % preallocates dy array
dy(:,:,1) = [-1 -1; 1 1 ]; 
dy(:,:,2) = [-1 -1; 1 1 ];
dy = 0.25 * dy;

dz = zeros(2,2,2); % preallocates dz array
dz(:,:,1) = [-1 -1; -1 -1 ]; 
dz(:,:,2) = ones(2,2);
dz = 0.25 * dz;

dt = ones(2,2,2); % creates a time differential
dt = 0.25*dt; 


[Ix,Iy,Iz] = gradient(image2,1,1,((8 * 13) /(55 * 1.25))); % takes the derivative with respect to time
It = image3-image1; 
It = It ./ (1 * (2 * delt));



end

