%==========================================================================
%% 3D Image Reader
%
% Read a series of images
%
% Original: Peshala T Gamage
% Refactored: Brendan Bouchard
% 20251125
% Last Updated: 20260112
%==========================================================================

clc
close all
clear 

s = 2;   %  slice of heart (2 to 11), select points for each slice changing s
n1 = 48 + (s * 30); % segment of slices to be used in for loop
n2 = n1 + 29; % unused

for p = n1:(n1 + 10)
    filename = sprintf('images/IM_%0.4d', p); % defines the file name for the image based on the iteration
    X = dicomread(filename); % reads the image file
    info = dicominfo(filename); % reads the metadata for the image file
    q = p + 1 - n1; % total number of slices to consider
    time(p - n1 + 1) = info.TriggerTime; % writes the time step to the metadata for the new image (maybe?)

    % yy=X(60:190,60:190);
    yy = X;
    % Y=imcrop(yy,[150 80 75 100]);
    % xx=imresize(Y,[541 541]);
    J(:,:,q) = imadjust(yy); % adjusts image intensity of the each slice
    name = ['frame' num2str(q) '.png']; % names the new image file
    
    imshow((J(:,:,q))) % displays the adjusted image in a figure
    Image = getframe(gcf); % captures axes as movie frame
    imwrite(J(:,:,q), name) % writes figure to new image file
end

J = im2single(J); % converts J to image with class single
