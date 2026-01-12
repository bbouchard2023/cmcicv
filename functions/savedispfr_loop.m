function [DX,DY,DZ]=savedispfr_loop(L,k,Points,start,stop,opt)
% Calculates displacement relative to the first frame

Ux=L.Ux;
Uy=L.Uy;
Uz=L.Uz;

DX{1}=zeros(1,length(Points)); % preallocates dx array
DY{1}=zeros(1,length(Points)); % preallocates dy array
DZ{1}=zeros(1,length(Points)); % preallocates dz array
dt=0.027/2; % 
t=dt;
% start=29; stop=59;  % start and end frams in loop (1:30 -> 1:30)
dx=zeros(length(Points),length(start:stop));
dy=zeros(length(Points),length(start:stop));
dz=zeros(length(Points),length(start:stop));
% fr_num=[1:30 1:30 1:30];
for loop=1:1 % does nothing
   
for jj=start:stop
%     fr=fr_num(jj);
    fr=jj % frames
   
    nn=jj-start+2 % 2 frames ahead of the current frame
%     Ux{1,fr}=imgaussfilt3(Ux{1,fr},0.25,'FilterSize',3);
%     Uy{1,fr}=imgaussfilt3(Uy{1,fr},0.25,'FilterSize',3);
%     Uz{1,fr}=imgaussfilt3(Uz{1,fr},3);
for i=1:length(Points)
    px(i)=Points{i}(1)-60;
    py(i)=Points{i}(2)-60;        % check whether -60 or -30  (where original image is cropped)
%     ceil(px(i)+dx(i,nn-1)*dt
    dx(i,nn)=dx(i,nn-1)+t*Ux{1,fr}(ceil(px(i)+dx(i,nn-1)*dt),ceil(py(i)+dy(i,nn-1)*dt),ceil(k+dz(i,nn-1)*dt));
    dy(i,nn)=dy(i,nn-1)+t*Uy{1,fr}(ceil(px(i)+dx(i,nn-1)*dt),ceil(py(i)+dy(i,nn-1)*dt),ceil(k+dz(i,nn-1)*dt));
    dz(i,nn)=dz(i,nn-1)+t*Uz{1,fr}(ceil(px(i)+dx(i,nn-1)*dt),ceil(py(i)+dy(i,nn-1)*dt),ceil(k+dz(i,nn-1)*dt));

    
end

DX(nn)={dx(:,nn)'};
DY(nn)={dy(:,nn)'};
DZ(nn)={dz(:,nn)'};


% scatter3(Points{i}(1),Points{i}(2),k)
% hold on
if opt==1 % not defined anywhere else
drawnow
scatter3(px+dy(:,nn)',py+dx(:,nn)',k+dz(:,nn)','filled','sc')
hold on
% plot3(px+dx,py+dy,k+dz)
axis([80 180 80 180 0 120]);
view([0 0 1]);


end
% scgx(nn-1)=(mean(dx(:,nn)));
% scgy(nn-1)=(mean(dy(:,nn)));
% scgz(nn-1)=(max(dz(:,nn)));
end
end
if opt==1
figure;
plot(scgx)
hold on
plot(scgy)
hold on
plot(scgz)
end