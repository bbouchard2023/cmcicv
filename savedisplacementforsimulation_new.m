%==========================================================================
%% Displacement for Simulation
%
% Saves displacements of tracked points
%
% Original: Peshala T Gamage
% Refactored: Brendan Bouchard
% 20251125
% Last Updated: 20251125
%==========================================================================

clc
clear
close all


vid = VideoWriter('Fina22.avi');
vid.FrameRate=60;
open(vid);
axis tight manual 
set(gca,'nextplot','replacechildren')

% load('time.mat')
L=load('New_1.mat');  % load tracked displacement from motion3Dfrbfr_tinterp_looped.m

clear X_disp  Y_disp  Z_disp
for loop=1:1
    
for t=1:61
    X=[]; Y=[]; Z=[]; eX=[]; eY=[]; eZ=[]; C=[]; x=[]; y=[]; z=[];
    Corx=[]; Cory=[]; Corz=[];
for s=3:10
    
    fname1=['EndoPoints' num2str(s) '.mat'];  % load poits to be tracked, for each slice
%     
    P1=load(fname1);

    Points1=P1.Points;
%     
    if s==1
    [DX,DY,DZ]=savedispfr_loop(L,s,Points1,1,60,0);

    else
     [DX,DY,DZ]=savedispfr_loop(L,round((s-1)*8/1.8909+1),Points1,1,60,0);
%      [DX,DY,DZ]=savedispfr_loop(L,round((s)),Points1,1,60,0);
%     
    end
    clear px epx
    clear py epy
    for i=1:length(Points1)
    px(i)=Points1{i}(1)-60;
    py(i)=Points1{i}(2)-60;
    end

    grid off
    plot3((px+DY{1,t}).*1.25,(py+DX{1,t}).*1.25,(s-1)*8+DZ{1,t}.*1.25,'k--','LineWidth',1.5)
    
%     
    hold on
    grid off
    xlim([(70-60)*1.25 (200-60)*1.25])
    ylim([(70-60)*1.25 (200-60)*1.25])
    zlim([5 100])
%     view([1 0 1])
    s
    t
    X=[(px+DY{1,t}).*1.25 X];
    Y=[(py+DX{1,t}).*1.25 Y];
    Z=[(s-1)*8+DZ{1,t}.*1.25 Z];
    x=[(DY{1,t}).*1.25 x];
    y=[(DX{1,t}).*1.25 y];
    z=[DZ{1,t}.*1.25 z];
    Corx=[(px).*1.25 Corx];
    Cory=[(py).*1.25 Cory];
    Corz=[s*8*ones(1,length(DZ{1,t})) Corz];
    
%     (px+60)'.*1.25 (py+60)'.*1.25 s*ones(length(px),1)];
end 
X_disp(t)={x'};
Y_disp(t)={y'};
Z_disp(t)={z'};

hold on
   cor3=[X(1:end);Y(1:end);Z(1:end)];

cor3d=cor3';


 DT = delaunayTriangulation(cor3d);  

%  tetramesh(DT,'FaceAlpha',0.3);
[K,v] = convexHull(DT);

V(t)=v;
trisurf(K,DT.Points(:,1),DT.Points(:,2),DT.Points(:,3),'EdgeColor','none','FaceColor','r','FaceAlpha',0.6,'LineWidth',0.01)
hold on
% trisurf(eK,eDT.Points(:,1),eDT.Points(:,2),eDT.Points(:,3),'EdgeColor','none','FaceColor','r','FaceAlpha',0.6,'LineWidth',0.01)
   xlim([(70-60)*1.25 (200-60)*1.25])
    ylim([(70-60)*1.25 (200-60)*1.25])
    zlim([5 80])
    view([1 0 1])
   frame = getframe(gcf);
   writeVideo(vid,frame);
   hold off
end
end
close(vid);

%% save displacement data at each point
for i=1:length(X_disp)
XDISP(:,i)=X_disp{1,i};
YDISP(:,i)=Y_disp{1,i};
ZDISP(:,i)=Z_disp{1,i};
end

CORD=[Corx' Cory' Corz'];

X_ansys=[CORD XDISP(:,2:end-1) XDISP(:,3:end)];
Y_ansys=[CORD YDISP(:,2:end-1) YDISP(:,3:end)];
Z_ansys=[CORD YDISP(:,2:end-1) YDISP(:,3:end)];

% x y z displacements for the simulation
X_ansys=X_ansys.*1.3;  
Y_ansys=Y_ansys.*1.3;
Z_ansys=Z_ansys.*1.3;