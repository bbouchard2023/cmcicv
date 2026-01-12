 %==========================================================================
%% Lucas-Kanade Optical Flow for Two 3D Images
%
% This function estimates deformations between two subsequent 3-D images
% using Lucas-Kanade optical flow equation. 
%
% Original: Peshala T Gamage
% Refactored: Brendan Bouchard
% 20251212
% Last Updated: 20251212
%==========================================================================

function [ux,uy,uz]=LK3D2_gauss_zvary( image1, image2,image3,r,d,fr,delt ) %   -image1, image2 :   two subsequent images or frames

%   -d : radius in z direction



 % ** in MATLAB the first dimension is y and 2nd is x and 3rd is y (Peshala)

%  Default parameter
if nargin==2
    r=2; %   -r : radius of the neighbourhood, default value is 2. 
end

[height,width,depth]=size(image1); 
image1=im2double(image1);
image2=im2double(image2);
image3=im2double(image3);
% Initializing flow vectors
ux=zeros(size(image1)); 
uy=ux; 
uz=ux;

% Computing image derivatives
[Ix,Iy,Iz,It]=imageDerivatives3D2(image1,image2,image3,fr,delt);
%% gaussian
% sigma=2.5;
% [x,y,z] = ndgrid(-r:r,-r:r,-r:r);
%     arg   = -(x.*x + y.*y + z.*z)/(2*sigma*sigma);
% w = exp(arg);
w=ones(length(-r:r),length(-r:r),length(-d:d)); % creates an identity matrix (kinda)
ww=w.*w; % not sure this actually does anything
%%
for i=(r+1):(height-r)
    for j=(r+1):(width-r)
        for k=(d+1):(depth-d)
        
        blockofIx=Ix(i-r:i+r,j-r:j+r,k-d:k+d); % takes surrounding values from point in Ix and puts them in an array
        blockofIy=Iy(i-r:i+r,j-r:j+r,k-d:k+d); % " in Iy
        blockofIz=Iz(i-r:i+r,j-r:j+r,k-d:k+d); % " in Iz
        blockofIt=It(i-r:i+r,j-r:j+r,k-d:k+d); % " in It

               
        A=zeros(3,3);
        B=zeros(3,1);
        
        A(1,1)=sum(sum(sum(blockofIx.^2.*ww))); % why three times?
        A(1,2)=sum(sum(sum(blockofIx.*blockofIy.*ww)));
        A(1,3)=sum(sum(sum(blockofIx.*blockofIz.*ww)));
        
        A(2,1)=sum(sum(sum(blockofIy.*blockofIx.*ww)));
        A(2,2)=sum(sum(sum(blockofIy.^2.*ww)));
        A(2,3)=sum(sum(sum(blockofIy.*blockofIz.*ww)));

        A(3,1)=sum(sum(sum(blockofIz.*blockofIx.*ww)));
        A(3,2)=sum(sum(sum(blockofIz.*blockofIy.*ww)));
        A(3,3)=sum(sum(sum(blockofIz.^2.*ww)));
       
        B(1,1)=sum(sum(sum(blockofIx.*blockofIt.*ww)));
        B(2,1)=sum(sum(sum(blockofIy.*blockofIt.*ww)));
        B(3,1)=sum(sum(sum(blockofIz.*blockofIt.*ww)));
        
        invofA=pinv(A); % pseudoinverse of A 
        
        V=invofA*(-B);
        ux(i,j,k)=V(1,1);
        uy(i,j,k)=V(2,1);
        uz(i,j,k)=V(3,1);
        end
    end
end

end

% Reference :
% Lucas, B. D., Kanade, T., 1981. An iterative image registration 
% technique with an application to stereo vision. In: Proceedings of the 
% 7th international joint conference on Artificial intelligence - Volume 2.
% Morgan Kaufmann Publishers Inc., San Francisco, CA, USA, pp. 674-679.
