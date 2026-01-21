clc
clear
close all

compmodel = readmatrix("SCG Displacement.csv");

t = compmodel(:,1);
displace = compmodel(:,2);

% [dispfreq,dispamp] = FFT_fn(disp,74.074);

%% Original Signal BPM and Peaks


[origpeak, xorig] = findpeaks(displace,"MinPeakHeight",0.4e-04);

orignumpeaks = length(origpeak);

deltatorig = zeros(1,length(origpeak));

for i = 2:length(xorig)
    deltatorig(i) = t(xorig(i)) - t(xorig(i - 1));
end
deltatorig(1) = [];
meandeltatorig = mean(deltatorig);
bpsorig = 1/meandeltatorig;
bpmorig = bpsorig * 60;

%% Filtered Signal BPM and Peaks
for j = 0.5:0.01:2
    for k = 2.01:0.01:5
        dispfilt = filter(filter_bandpass_disp(j,k),displace);

        [peak, x] = findpeaks(dispfilt,"MinPeakHeight",0.2e-04);

        numpeaks = length(peak);

        deltat = zeros(1,length(peak));

        for i = 2:length(x)
            deltat(i) = t(x(i)) - t(x(i - 1));
        end
        meandeltat = mean(deltat);
        bps = 1/meandeltat;
        bpm = bps * 60;

        deltapeaks = abs(numpeaks - orignumpeaks);
        deltabpm = abs(bpm - bpmorig);

        fprintf("j: %d k: %d\nDifference in BPM: %f\n\n",j,k,deltabpm)
        if deltabpm < 20 && deltabpm > 0
            break
        else
            continue
        end
    end
end


%% Plotting

fig = figure;
fig.WindowState = 'maximize';
plot(t,displace,"b")
hold on
plot(t,dispfilt,"g")
scatter(t(x),peak)
scatter(t(xorig),origpeak,"m")
title("ANSYS Displacement")
xlabel("Time (s)")
ylabel("Amplitude")
legend("Raw Displacement", "Bandpass Filtered Displacement","Peaks","Unfiltered Peaks")
hold off

fprintf("For filtered:\nNumber of Peaks: %d\nBPM: %f\n\n",numpeaks,bpm)
fprintf("For original:\nNumber of Peaks: %d\nBPM: %f\n\n",orignumpeaks,bpmorig)
fprintf("Difference in Peaks: %d\nDifference in BPM: %f\n",deltapeaks,deltabpm)

