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

%% This is going to need a magician to decode

function [Ix, Iy, Iz, It] = imageDerivatives3D2(image1, image2, image3, fr, delt)


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

dt = ones(2,2,2);
dt = 0.25*dt;

% Computing derivatives
% Ix = 0.5 * (convn(image1,dx) + convn(image2,dx) );
% Iy = 0.5 * (convn(image1,dy) + convn(image2,dy) );
% Iz = 0.5 * (convn(image1,dz) + convn(image2,dz) );
% It = 0.5* (convn(image2,dt) - convn(image1,dt) );
% [Ix,Iy,Iz]=gradient(image1,1,1,8);
% if fr==2 || fr==30
%     It=image2-image1;
%     It=It./delt;
% end
% Ix = 1* (convn(image2,dx));
% Iy = 1 * (convn(image2,dy));
% Iz = (1/6) * (convn(image2,dz));
% It = (1/2*delt)* (convn(image3,dt) - convn(image1,dt) );  
[Ix,Iy,Iz] = gradient(image2,1,1,((8 * 13) /(55 * 1.25)));
It = image3-image1;
It = It ./ (1 * (2 * delt));
% It=image2-image1;
% Adjusting sizes
% Ix=Ix(1:size(Ix,1)-1, 1:size(Ix,2)-1, 1:size(Ix,3)-1);
% Iy=Iy(1:size(Iy,1)-1, 1:size(Iy,2)-1, 1:size(Iy,3)-1);
% Iz=Iz(1:size(Iz,1)-1, 1:size(Iz,2)-1, 1:size(Iz,3)-1);
% It=It(1:size(It,1)-1, 1:size(It,2)-1, 1:size(It,3)-1);


end

