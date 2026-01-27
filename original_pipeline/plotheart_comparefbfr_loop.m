%==========================================================================
%% 3D Motion Visualizer
%
% Visualizes movement of tracked points on image  
%
% Original: Peshala T Gamage
% Refactored: Brendan Bouchard
% 20251125
% Last Updated: 20260112
%==========================================================================

clc
clear 


clear Points X_mov Y_mov
load('Timeinterpdata.mat') % image data
load('New_1.mat')  % load tracked displacement

load('HPoints9.mat') % load points
k=9;  % slice number (should match with the selected points)

n1=48+k*30; % unused
clear dx px
clear dy py
clear dz 
clear dz 
clear scgx(fr-1) scgy(fr-1) scgz(fr-1)

vid = VideoWriter('111.avi');  % video name
vid.FrameRate=30; % frame rate
open(vid); % opens video
axis tight manual 
set(gca,'nextplot','replacechildren')

DX{1}=zeros(1,length(Points)); % preallocates array for dx
DY{1}=zeros(1,length(Points)); % " for dy
DZ{1}=zeros(1,length(Points)); % " for dz
dt=0.027/2; % creates time differential
t=dt;
start=2; stop=58;  % start and end frames in loop (1:30 -> 1:30)
dx=zeros(length(Points),length(start:stop));
dy=zeros(length(Points),length(start:stop));
dz=zeros(length(Points),length(start:stop));
% fr_num=[1:30 1:30 1:30];

for loop=1:1
   
for jj=start:stop % for frames start - stop
%     fr=fr_num(jj);
    fr=jj; % defines frames as fr and jj
   
    nn=jj-start+2 % 2 frames ahead of the current frame
    %% gauss smoothing to velocity field
    Ux{1,fr}=imgaussfilt3(Ux{1,fr},0.25,'FilterSize',3);
    Uy{1,fr}=imgaussfilt3(Uy{1,fr},0.25,'FilterSize',3);
    Uz{1,fr}=imgaussfilt3(Uz{1,fr},0.5,'FilterSize',3);
    %% sets the dx, dy, and dz for each selected point
for i=1:length(Points)
    px(i)=Points{i}(1)-60;
    py(i)=Points{i}(2)-60;
    dx(i,nn)=dx(i,nn-1)+t*Ux{1,fr}(ceil(px(i)+dx(i,nn-1)*dt),ceil(py(i)+dy(i,nn-1)*dt),ceil(k*8/1.8909+1+dz(i,nn-1)*dt));
    dy(i,nn)=dy(i,nn-1)+t*Uy{1,fr}(ceil(px(i)+dx(i,nn-1)*dt),ceil(py(i)+dy(i,nn-1)*dt),ceil(k*8/1.8909+1+dz(i,nn-1)*dt));
    dz(i,nn)=dz(i,nn-1)+t*Uz{1,fr}(ceil(px(i)+dx(i,nn-1)*dt),ceil(py(i)+dy(i,nn-1)*dt),ceil(k*8/1.8909+1+dz(i,nn-1)*dt));

  
    
end

DX(nn)={dx(:,nn)'}; % sets dx in the DX cell array
DY(nn)={dy(:,nn)'}; % sets dy in the DY cell array
DZ(nn)={dz(:,nn)'}; % sets dz in the DZ cell array

% scatter3(Points{i}(1),Points{i}(2),k)
hold off
drawnow
scatter3(px+dy(:,nn)'.*1.25,py+dx(:,nn)'.*1.25,8+dz(:,nn)'.*1.15,100,'filled','sc')
hold on
plot3(px+dy(:,nn)',py+dx(:,nn)',8+dz(:,nn)','y','LineWidth',2.5)
axis([80 180 80 180 0 120]);
hold on
quiver3(px+dy(:,nn)'.*1.25,py+dx(:,nn)'.*1.25,8+dz(:,nn)'.*1.15,dy(:,nn)'.*1,dx(:,nn)'.*1,dz(:,nn)'.*0, 0.65,'r','LineWidth',2.5) % vectors
hold on
% scatter3(X_mov(:,1)-60+X_mov(:,fr+1),Y_mov(:,2)-60+Y_mov(:,fr+1),k*ones(1,length(Points)),'filled','r')

axis([80-60 180-60 80-60 180-60 0 100]);
view([0 0 1]);

set(gcf,'Position',[100 100 500 500])
max(dz);

%% plotting the images
% X=Iinterp(:,:,round(k*8/1.8909+1),fr); 
X=Iinterp(:,:,round(k),fr); 

yy=X(60:end,60:end);   % crop image
zzz=yy;
hold on
imshow(imadjust(zzz))  % show image

view([0 0 1])
set(gcf,'Position',[100 100 500 500])

%% Rigid body displacement

scgx(nn-1)=(mean(dx(:,nn)));  % mean displacement x
scgy(nn-1)=(mean(dy(:,nn))); % mean displacement y
scgz(nn-1)=(max(dz(:,nn)));% mean displacement z

set(gcf,'Position',[100 100 500 500])
frame = getframe(gcf);
   writeVideo(vid,frame);
   hold off
end
end
close(vid);
figure;
plot(scgx)
hold on
plot(scgy)
hold on
plot(scgz)