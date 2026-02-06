%==========================================================================
%% Displacement Comparison
%
% Compares displacement at heart wall and chest
%
% Brendan Bouchard
% 20260206
% Last Updated: 20260206
%==========================================================================

clc
clear
close all

dispdata = readmatrix("External_Internal_Displacement.csv");
t_ext = dispdata(:,1);
disp_ext = dispdata(:,2);

t_int = dispdata(:,4);
disp_int = dispdata(:,5);

intzeros = [];
extzeros = [];
diffint = diff(disp_int);
diffext = diff(disp_ext);

for i = 2:length(diffint)
    if diffint(i) < 0 && diffint(i - 1) > 0
        intzeros = cat(1,intzeros,[t_int(i) disp_int(i)]);
    end
end

for i = 2:length(diffext)
    if diffext(i) < 0 && diffext(i - 1) > 0
        extzeros = cat(1,extzeros,[t_ext(i) disp_ext(i)]);
    end
end

displag = abs(extzeros(:,1) - intzeros(:,1));


subplot(2,1,1)
yyaxis left
plot(t_int,disp_int,"m")
hold on
scatter(intzeros(:,1),intzeros(:,2),"yx")
yyaxis right
plot(t_ext,disp_ext,"g")
scatter(extzeros(:,1),extzeros(:,2),"b")
xlabel("Time (s)")
ylabel("Displacement (mm)")
legend("Internal Displacement","Internal Peaks","External Displacement","External Peaks")
title("Displacement Comparison")
hold off

subplot(2,1,2)
plot(displag,"r")
xlabel("Displacement Peak")
ylabel("Difference in Peak Time")
title("Lag Between Cardiac Wall Displacement and Chest Displacement")

saveas(gcf,"DisplacementComparison.png")
