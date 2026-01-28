clear all
% to interpolate the 4D data in time

for i=1:30
 
Vname=['volumeinterp' num2str(i) '.mat'];
load(Vname)
I(:,:,:,i)=J_new  ;

end
I(:,:,:,31)=I(:,:,:,1);

F = griddedInterpolant(I,'spline');
x=(1:1:256);
y=(1:1:256);
z=(1:1:           56);
tq=(1:0.5:31);

Iinterp=F({x,y,z,tq});

save('Timeinterpdata_2.mat','Iinterp')