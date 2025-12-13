%==========================================================================
%% Point Selector
%
% This code is used to select points, points should be selected for
% Each slice s is selected in dicom_read
%
% Original: Peshala T Gamage
% Refactored: Brendan Bouchard
% 20251125
% Last Updated: 20251125
%==========================================================================


clc
clear
close all

% X1=K;
k = 1;          % frame number
I = J(:,:,k);   % image stack
% I(0.5<I<0.95)=0;
% I = imgaussfilt(I,2);
I = cat(3,I,I,I);  % to make the image compatiable with the code (R G B)
clear xi
clear yi

figure;  
imshow(I); % displays image in figure
title('Select Points');  % pick the points using mouse clicks, press right click for the last point
[xi,yi] = getpts; % specify points by mouse click (for boundary selection)
Pname = ['EndoPoints_new' num2str(s) '.mat'];  % give a name 
% load(Pname)
pp = []; % creates an array to store the points for each selection
Points = {}; % creates a cell array to save the points
for i = 1:length(xi)
    pp(i,1) = xi(i); % saves the X component of each point
    pp(i,2) = yi(i); % saves the Y component of each point
    Points{i} = pp(i,:); % writes the array into a cell array
end
hold on
a = (1:length(xi))'; 
b = num2str(a); 
c = cellstr(b);

smx = 0.5; % displacement so the text does not overlay the data points
smy = 0.1; % displacement so the text does not overlay the data points
scatter(xi,yi,'filled') % scatter plot of something or other
text(xi+smx, yi+smy, c,'Color','y'); 

save(Pname,'Points') % saves Points cell array to file