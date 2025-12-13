%==========================================================================
%% 3D Motion Calculator
%
% Calculates 3D motion, takes lot of time to run. use parallel processing
% toolbox if possible, then change the for loops to 'parfor'  
%
% Original: Peshala T Gamage
% Refactored: Brendan Bouchard
% 20251125
% Last Updated: 20251125
%==========================================================================


clc
clear 

load('Timeinterpdata_2.mat');  % 4D image data ( 3D data in time), see the folder original_data to see how data is prepared


for loop=1:2
    for fr=1:60 % frames 1-60
        fr_l=fr+(loop-1)*60; % frame rate correction
        if fr==1
            img1=Iinterp(60:180,30:180,:,60); % defines img1 as a range inside the array Iinterp1 (from .mat file)
        else
            img1=Iinterp(60:180,30:180,:,fr-1); % defines img1, the image before the one in iteration
        end
        img2=Iinterp(60:180,30:180,:,fr); % defines the image in the iteration
        img3=Iinterp(60:180,30:180,:,fr+1); % defines the image after the one in iteration


        % img1 = imgaussfilt3(img1);
        % img2 = imgaussfilt3(img2);
        % img3 = imgaussfilt3(img3);
        for i=1:size(img1,3)
            % adjusts intensity of pixels for use in motion tracking (probably increases contrast)
            img1(:,:,i)=imadjust(img1(:,:,i)); 
            % img1(:,:,i)=imdiffusefilt(img1(:,:,i));
            img2(:,:,i)=imadjust(img2(:,:,i));
            % img2(:,:,i)=imdiffusefilt(img2(:,:,i));
            img3(:,:,i)=imadjust(img3(:,:,i));
            % img3(:,:,i)=imdiffusefilt(img3(:,:,i));
        end

        [ux,uy,uz]=LK3D2_gauss_zvary( img1, img2,img3,10,5,fr,0.027/2);  % motion tracking function


        Ux(fr_l)={ux}; % defines where to save motion tracking data in X direction
        Uy(fr_l)={uy}; % " in Y direction
        Uz(fr_l)={uz}; % " in Z direction



    end
end
save('New_1.mat' ,'Ux','Uy','Uz');  % save the displacement field





