%% ============================================================
%  SCG/Displacement Respiration Bridge Pipeline (ANSYS -> MATLAB)
%  Goal: Extract respiration from displacement and show it is
%  recoverable from SCG-like acceleration via physics-consistent
%  integration/differentiation + robust filtering.
%  Author: Generated Script   Date: January 19, 2026
%
%  KEY IMPROVEMENTS:
%  1) Uses cumtrapz(t, signal) - proper time-vector integration
%  2) Sanity check: differentiates displacement to verify consistency
%  3) Two integration methods: polynomial detrending & FFT-based
%  4) Separates BRIDGE (acc↔disp) from FEATURE extraction (cardiac/resp)
%  5) Paper-friendly metrics: cardiac waveform + resp rate validation
%% ============================================================

clear; clc; close all;

%% ---------------------- 0) CONFIG ---------------------------
CFG = struct();

% Files (edit)
CFG.data.dispFile = "SCG Displacement.csv";     % time, disp (m) or (mm)
CFG.data.accFile  = "SCG Acceleration.csv";     % time, acc  (m/s^2)
CFG.data.hasSeparateTimeColumns = true;

% Sampling
CFG.fs = [];                        % leave [] to estimate from time vector

% Units
CFG.units.disp = "m";               % "m" or "mm"
CFG.units.acc  = "m/s^2";

% Respiration band (adult typical)
CFG.resp.f_lo = 0.10;               % Hz
CFG.resp.f_hi = 0.50;               % Hz

% Windowing for resp rate tracking
CFG.resp.winSec = 30;               % seconds (20–40 is common)
CFG.resp.hopSec = 5;                % seconds

% SCG-like band (for morphology/envelope) - adjust as needed
CFG.scg.f_lo = 5;                   % Hz
CFG.scg.f_hi = 35;                  % Hz (must be < Nyquist frequency)

% Drift control for integration (very low cutoff)
CFG.intDriftHP = 0.05;              % Hz (must be < resp.f_lo to preserve resp)

% Differentiation noise control
CFG.diff.sgolayOrder = 3;
CFG.diff.sgolayFrameSec = 0.25;     % seconds (0.15–0.35 typical)

% Plot toggles
CFG.plot.makeFigures = true;

%% ---------------------- 1) LOAD DATA ------------------------
fprintf('Loading ANSYS data...\n');
D = load_ansys_signals(CFG);

% Update CFG.fs if it was estimated
if isempty(CFG.fs)
    CFG.fs = D.fs_est;
end

% Optional: sanity checks
assert(all(isfinite(D.t)) && all(isfinite(D.disp)) && all(isfinite(D.acc)), ...
    "Non-finite values found in input signals.");

fprintf('Data loaded: %d samples, fs = %.2f Hz, duration = %.2f s\n', ...
    length(D.t), CFG.fs, D.t(end)-D.t(1));

% Adjust window size if signal is too short
signal_duration = D.t(end) - D.t(1);
if signal_duration < CFG.resp.winSec
    fprintf('WARNING: Signal duration (%.2f s) is shorter than respiration window (%.2f s)\n', ...
        signal_duration, CFG.resp.winSec);
    CFG.resp.winSec = max(2, signal_duration * 0.8);  % Use 80% of signal length
    CFG.resp.hopSec = max(0.5, CFG.resp.winSec * 0.2); % 20% overlap
    fprintf('Adjusted window to %.2f s, hop to %.2f s\n', CFG.resp.winSec, CFG.resp.hopSec);
end

% Check and adjust SCG band to not exceed Nyquist frequency
nyquist_freq = CFG.fs / 2;
if CFG.scg.f_hi >= nyquist_freq
    CFG.scg.f_hi = nyquist_freq * 0.95;  % Set to 95% of Nyquist
    fprintf('WARNING: SCG upper frequency adjusted to %.2f Hz (below Nyquist = %.2f Hz)\n', ...
        CFG.scg.f_hi, nyquist_freq);
end

% Note about signal duration for respiration analysis
if signal_duration < 20
    fprintf('NOTE: Signal duration (%.2f s) is very short for respiration analysis.\n', signal_duration);
    fprintf('      For reliable respiration rate tracking, signals of 30-60+ seconds are recommended.\n');
end

%% ---------------------- 2) PREPROCESS -----------------------
fprintf('Preprocessing signals...\n');
% 2a) Unit normalization
if CFG.units.disp == "mm"
    D.disp = D.disp * 1e-3;
end

% 2b) Remove mean (good practice before integration)
D.disp0 = detrend(D.disp, 'constant');
D.acc0  = detrend(D.acc,  'constant');

% 2c) Optional: de-spike / clip extreme outliers
% D.disp0 = filloutliers(D.disp0,'linear','movmedian',round(CFG.fs*0.5));
% D.acc0  = filloutliers(D.acc0 ,'linear','movmedian',round(CFG.fs*0.5));

%% ----------- 3) SANITY CHECK: DISP -> ACC (DIFFERENTIATION) -----
fprintf('Sanity check: differentiating displacement twice...\n');
% This verifies that acceleration is actually the second derivative of displacement
D.acc_from_disp = disp_to_acc(D.disp0, D.t, CFG.diff);

% Compare with actual acceleration
acc_corr = corrcoef(D.acc0, D.acc_from_disp);
D.sanity_check_corr = acc_corr(1,2);
fprintf('  Correlation (actual acc vs diff of disp): %.4f\n', D.sanity_check_corr);
if D.sanity_check_corr < 0.7
    fprintf('  WARNING: Low correlation suggests signals may not be from same node/axis!\n');
end

%% ----------- 4) BRIDGE: ACC -> DISP (DOUBLE INTEGRATION) -----
fprintf('Computing displacement from acceleration (double integration)...\n');
% Key: integrate with drift control WITHOUT killing cardiac/respiration bands.

% Option A: Time-domain integration with polynomial drift removal
fprintf('  Option A: Polynomial detrending...\n');
D.disp_from_acc_polyfit = acc_to_disp_polyfit(D.acc0, D.t);

% Option B: FFT-based integration (avoids time-domain drift)
fprintf('  Option B: FFT-based integration...\n');
D.disp_from_acc_fft = acc_to_disp_fft(D.acc0, D.t, CFG.fs, 0.01);

% Choose method based on which produces valid results
fprintf('  Polyfit result: range=[%.2e, %.2e], finite=%d/%d\n', ...
    min(D.disp_from_acc_polyfit), max(D.disp_from_acc_polyfit), ...
    sum(isfinite(D.disp_from_acc_polyfit)), length(D.disp_from_acc_polyfit));
fprintf('  FFT result:     range=[%.2e, %.2e], finite=%d/%d\n', ...
    min(D.disp_from_acc_fft), max(D.disp_from_acc_fft), ...
    sum(isfinite(D.disp_from_acc_fft)), length(D.disp_from_acc_fft));

if all(isfinite(D.disp_from_acc_fft))
    D.disp_from_acc = D.disp_from_acc_fft;
    fprintf('  Using FFT method for reconstruction\n');
else
    D.disp_from_acc = D.disp_from_acc_polyfit;
    fprintf('  FFT method failed, using polynomial detrending method\n');
end

% Final sanity check
if ~all(isfinite(D.disp_from_acc))
    error(['Integration failed to produce finite values. This usually means:\n' ...
           '  1) Input acceleration has non-finite values\n' ...
           '  2) Signals are from different nodes/axes (corr=%.4f < 0.7)\n' ...
           '  3) Acceleration includes rigid-body motion not in displacement\n' ...
           '  Please verify your ANSYS export settings.'], D.sanity_check_corr);
end

%% ---------------------- 5) CARDIAC & RESP EXTRACTION -------------------
fprintf('Extracting cardiac and respiration components...\n');

% CARDIAC BAND (5-40 Hz) - for waveform similarity validation
fprintf('  Cardiac band (5-40 Hz)...\n');
C.disp_cardiac = bandpass_zerophase(D.disp0, CFG.fs, CFG.scg.f_lo, CFG.scg.f_hi);
C.dispA_cardiac = bandpass_zerophase(D.disp_from_acc, CFG.fs, CFG.scg.f_lo, CFG.scg.f_hi);

% RESPIRATION BAND (0.1-0.5 Hz) - for respiration rate tracking
fprintf('  Respiration band (0.1-0.5 Hz)...\n');
% 5a) Resp from displacement (direct)
R.disp_resp = extract_resp_component(D.disp0, CFG.fs, CFG.resp);
R.disp_rr   = track_resp_rate(R.disp_resp, CFG.fs, CFG.resp);

% 5b) Resp from reconstructed displacement (acc -> disp)
R.dispA_resp = extract_resp_component(D.disp_from_acc, CFG.fs, CFG.resp);
R.dispA_rr   = track_resp_rate(R.dispA_resp, CFG.fs, CFG.resp);

% 5c) Resp from SCG-like acceleration via envelope
A_scg = bandpass_zerophase(D.acc0, CFG.fs, CFG.scg.f_lo, CFG.scg.f_hi);
A_env = abs(hilbert(A_scg));
R.acc_env_resp = extract_resp_component(A_env, CFG.fs, CFG.resp);
R.acc_env_rr   = track_resp_rate(R.acc_env_resp, CFG.fs, CFG.resp);

%% ---------------------- 6) METRICS ---------------------------
fprintf('Computing metrics...\n');
M = struct();

% === PAPER-FRIENDLY VALIDATION METRICS ===

% 1) CARDIAC-BAND WAVEFORM SIMILARITY (Bridge validation)
fprintf('  1) Cardiac waveform similarity (5-40 Hz)...\n');
[M.cardiac_corr, M.cardiac_lagSec] = max_xcorr(C.disp_cardiac, C.dispA_cardiac, CFG.fs, 1);
M.cardiac_rmse = rmse_metric(C.disp_cardiac, C.dispA_cardiac);
M.cardiac_nrmse = M.cardiac_rmse / (max(C.disp_cardiac) - min(C.disp_cardiac)) * 100; % Normalized RMSE (%)

% 2) RESPIRATION-BAND AGREEMENT (Feature extraction validation)
fprintf('  2) Respiration component agreement (0.1-0.5 Hz)...\n');
[M.resp_corr, M.resp_lagSec] = max_xcorr(R.disp_resp, R.dispA_resp, CFG.fs, 10);

% 3) RESPIRATION RATE METRICS (Clinical metric validation)
fprintf('  3) Respiration rate tracking error...\n');
[M.rr_mae_bpm, M.rr_rmse_bpm] = rr_error_metrics(R.disp_rr.bpm, R.dispA_rr.bpm);

% 4) SCG-ENVELOPE vs DISPLACEMENT RR (Optional: alternative method)
[M.env_rr_mae_bpm, M.env_rr_rmse_bpm] = rr_error_metrics(R.disp_rr.bpm, R.acc_env_rr.bpm);

% 5) BROADBAND DISPLACEMENT MATCH (Overall reconstruction quality)
[M.broadband_corr, M.broadband_lagSec] = max_xcorr(D.disp0, D.disp_from_acc, CFG.fs, 1);
M.broadband_rmse = rmse_metric(D.disp0, D.disp_from_acc);
M.broadband_nrmse = M.broadband_rmse / (max(D.disp0) - min(D.disp0)) * 100;

fprintf('\n========================================\n');
fprintf('     VALIDATION METRICS SUMMARY\n');
fprintf('========================================\n\n');

fprintf('=== 1) CARDIAC BAND (5-40 Hz) Bridge Validation ===\n');
fprintf('Correlation:   %.4f\n', M.cardiac_corr);
fprintf('Lag:           %.4f seconds\n', M.cardiac_lagSec);
fprintf('NRMSE:         %.2f%%\n\n', M.cardiac_nrmse);

fprintf('=== 2) RESPIRATION BAND (0.1-0.5 Hz) Agreement ===\n');
fprintf('Correlation:   %.4f\n', M.resp_corr);
fprintf('Lag:           %.4f seconds\n\n', M.resp_lagSec);

fprintf('=== 3) RESPIRATION RATE Metrics ===\n');
fprintf('Disp RR vs Reconstructed Disp RR:\n');
fprintf('  MAE:  %.2f breaths/min\n', M.rr_mae_bpm);
fprintf('  RMSE: %.2f breaths/min\n\n', M.rr_rmse_bpm);

fprintf('Disp RR vs SCG-Envelope RR:\n');
fprintf('  MAE:  %.2f breaths/min\n', M.env_rr_mae_bpm);
fprintf('  RMSE: %.2f breaths/min\n\n', M.env_rr_rmse_bpm);

fprintf('=== 4) BROADBAND DISPLACEMENT Match ===\n');
fprintf('Correlation:   %.4f\n', M.broadband_corr);
fprintf('Lag:           %.4f seconds\n', M.broadband_lagSec);
fprintf('NRMSE:         %.2f%%\n\n', M.broadband_nrmse);

fprintf('========================================\n\n');

%% ---------------------- 7) PLOTS -----------------------------
if CFG.plot.makeFigures
    fprintf('Generating plots...\n');
    make_summary_plots(D, R, C, CFG, M);
end

%% ---------------------- 8) EXPORT ----------------------------
fprintf('Exporting results...\n');
% Save key results for paper figures/tables
OUT = struct("CFG",CFG,"D",D,"R",R,"C",C,"M",M);
save("scg_disp_resp_pipeline_output.mat","-struct","OUT");

% Export metrics to CSV
metrics_table = struct2table(M);
writetable(metrics_table,"metrics_summary.csv");

fprintf('Results saved to:\n');
fprintf('  - scg_disp_resp_pipeline_output.mat\n');
fprintf('  - metrics_summary.csv\n\n');

fprintf('========================================\n');
fprintf('     INTERPRETATION GUIDE\n');
fprintf('========================================\n\n');
fprintf('GOOD RECONSTRUCTION if:\n');
fprintf('  - Cardiac corr > 0.85 (waveform match)\n');
fprintf('  - Broadband corr > 0.80 (overall match)\n');
fprintf('  - Resp RR MAE < 2 bpm (rate accuracy)\n');
fprintf('  - NRMSE < 10%% (normalized error)\n\n');
fprintf('If metrics are poor, check:\n');
fprintf('  - Same node/axis in ANSYS export?\n');
fprintf('  - Acceleration type (absolute vs relative)?\n');
fprintf('  - Coordinate frame consistency?\n');
fprintf('  - Signal duration (need 30-60s for resp)\n');
fprintf('========================================\n\n');

disp('Pipeline complete!');

%% ============================================================
%                    LOCAL FUNCTIONS
%% ============================================================

function D = load_ansys_signals(CFG)
    % Expect CSV with columns: time, signal
    dispTab = readtable(CFG.data.dispFile);
    accTab  = readtable(CFG.data.accFile);

    t_disp = dispTab{:,1};
    x_disp = dispTab{:,2};

    t_acc  = accTab{:,1};
    x_acc  = accTab{:,2};

    % Align by interpolation onto a common time base
    t0 = max([t_disp(1), t_acc(1)]);
    t1 = min([t_disp(end), t_acc(end)]);

    % Choose denser time base (or use disp time)
    t = t_disp(t_disp>=t0 & t_disp<=t1);

    x_disp = interp1(t_disp, x_disp, t, 'linear', 'extrap');
    x_acc  = interp1(t_acc , x_acc , t, 'linear', 'extrap');

    % Estimate fs if not provided
    dt = median(diff(t));
    fs = 1/dt;

    D = struct();
    D.t = t(:);
    D.fs_est = fs;
    D.disp = x_disp(:);
    D.acc  = x_acc(:);
end

function y = bandpass_zerophase(x, fs, f1, f2)
    % Check for finite input
    if ~all(isfinite(x))
        error('Expected input signal to be finite. Found %d non-finite values.', sum(~isfinite(x)));
    end
    
    d = designfilt('bandpassiir','FilterOrder',4, ...
        'HalfPowerFrequency1',f1,'HalfPowerFrequency2',f2, ...
        'SampleRate',fs);
    y = filtfilt(d, x);
end

function y = highpass_zerophase(x, fs, fc)
    d = designfilt('highpassiir','FilterOrder',4, ...
        'HalfPowerFrequency',fc,'SampleRate',fs);
    y = filtfilt(d, x);
end

function disp_rec = acc_to_disp_polyfit(acc, t)
    % Option A: Time-domain integration with polynomial drift removal
    % Uses actual time vector (not /fs) and removes drift with polyfit
    
    % Ensure column vectors
    acc = acc(:);
    t = t(:);
    
    % Remove mean first
    acc = acc - mean(acc);
    
    % First integration: acc -> velocity
    v = cumtrapz(t, acc);
    
    % Remove polynomial drift from velocity (order 2)
    p = polyfit(t, v, 2);
    v_trend = polyval(p, t);
    v = v - v_trend;
    
    % Second integration: velocity -> displacement
    disp_rec = cumtrapz(t, v);
    
    % Remove polynomial drift from displacement (order 2)
    p = polyfit(t, disp_rec, 2);
    d_trend = polyval(p, t);
    disp_rec = disp_rec - d_trend;
    
    % Final mean removal
    disp_rec = disp_rec - mean(disp_rec);
end

function disp_rec = acc_to_disp_fft(acc, t, fs, hp_cutoff)
    % Option B: FFT-based integration (avoids time-domain drift)
    % Integrates in frequency domain with DC protection
    
    N = length(acc);
    
    % Remove mean first
    acc = acc - mean(acc);
    
    % FFT of acceleration
    ACC = fft(acc);
    
    % Frequency vector (proper fftshift-style)
    if mod(N, 2) == 0
        f = [(0:N/2), (-N/2+1:-1)] * (fs/N);
    else
        f = [(0:(N-1)/2), (-(N-1)/2:-1)] * (fs/N);
    end
    
    % Initialize displacement FFT
    DISP = zeros(size(ACC));
    
    % Integrate all frequencies except DC and very low frequencies
    for k = 1:N
        if abs(f(k)) >= hp_cutoff
            % Integrate twice: divide by (j*2*pi*f)^2 = -(2*pi*f)^2
            omega = 2*pi*f(k);
            DISP(k) = ACC(k) / (-omega^2);
        else
            % Kill DC and very low frequencies
            DISP(k) = 0;
        end
    end
    
    % Inverse FFT
    disp_rec = real(ifft(DISP));
    
    % Remove any remaining mean and ensure finite
    disp_rec = disp_rec - mean(disp_rec);
    
    % Sanity check
    if any(~isfinite(disp_rec))
        warning('FFT integration produced non-finite values. Falling back to polynomial method.');
        disp_rec = nan(size(acc));
    end
end

function acc_rec = disp_to_acc(dispSig, t, diffCFG)
    % Smooth displacement then differentiate twice using actual time vector
    % Uses Savitzky-Golay filter to smooth before differentiation
    
    fs = 1 / mean(diff(t));
    frame = max(5, round(diffCFG.sgolayFrameSec*fs));
    if mod(frame,2)==0, frame = frame+1; end % frame must be odd
    disp_s = sgolayfilt(dispSig, diffCFG.sgolayOrder, frame);

    % Use actual time vector for differentiation (not *fs)
    v = gradient(disp_s) ./ gradient(t);      % 1st derivative
    a = gradient(v) ./ gradient(t);           % 2nd derivative
    acc_rec = a;
end

function resp = extract_resp_component(x, fs, respCFG)
    % Isolate respiration band with zero-phase filtering
    d = designfilt('bandpassiir','FilterOrder',4, ...
        'HalfPowerFrequency1',respCFG.f_lo, ...
        'HalfPowerFrequency2',respCFG.f_hi, ...
        'SampleRate',fs);
    resp = filtfilt(d, x);
end

function RR = track_resp_rate(respSig, fs, respCFG)
    % Sliding-window Welch peak tracking in resp band
    win = round(respCFG.winSec * fs);
    hop = round(respCFG.hopSec * fs);
    n = length(respSig);
    
    % Check if signal is long enough for at least one window
    if n < win
        fprintf('WARNING: Signal too short for respiration rate tracking window.\n');
        fprintf('         Using entire signal as one window.\n');
        win = n;
        hop = n; % Only one window
    end

    idxStart = 1:hop:(n-win+1);
    
    % Ensure at least one window
    if isempty(idxStart)
        idxStart = 1;
    end
    bpm = nan(size(idxStart));
    tmid = nan(size(idxStart));

    for k = 1:numel(idxStart)
        i1 = idxStart(k);
        i2 = min(i1 + win - 1, n);

        xw = respSig(i1:i2);

        % Welch PSD - adjust window size for short signals
        welch_win = min(round(length(xw)*0.5), length(xw)-1);
        welch_win = max(welch_win, 8); % Minimum window size
        [Pxx,F] = pwelch(xw, hamming(welch_win), [], [], fs);

        % Find peak in resp band
        bandMask = (F>=respCFG.f_lo) & (F<=respCFG.f_hi);
        Fb = F(bandMask); Pb = Pxx(bandMask);

        [~,im] = max(Pb);
        fpeak = Fb(im);

        bpm(k) = fpeak * 60;
        tmid(k) = (i1+i2)/2 / fs;
    end

    RR = struct("tmid",tmid(:),"bpm",bpm(:));
end

function [cmax, lagSec] = max_xcorr(x, y, fs, maxLagSec)
    % Max normalized cross-correlation within +/- maxLagSec
    maxLag = round(maxLagSec*fs);
    [c,lags] = xcorr(zscore(x), zscore(y), maxLag, 'coeff');
    [cmax, im] = max(c);
    lagSec = lags(im)/fs;
end

function [mae, rmse] = rr_error_metrics(bpm1, bpm2)
    % Compare two BPM tracks (same length assumed; if not, align before)
    n = min(numel(bpm1), numel(bpm2));
    e = bpm1(1:n) - bpm2(1:n);
    mae = mean(abs(e),'omitnan');
    rmse = sqrt(mean(e.^2,'omitnan'));
end

function rmse_val = rmse_metric(x, y)
    % Root mean square error between two signals
    n = min(length(x), length(y));
    e = x(1:n) - y(1:n);
    rmse_val = sqrt(mean(e.^2));
end

function make_summary_plots(D, R, C, CFG, M)
    t = D.t;
    fs = CFG.fs;

    % Figure 1: Broadband signals comparison (Bridge validation)
    figure('Name', 'Bridge Validation - Broadband Signals', 'Position', [100, 100, 1400, 900]);
    
    subplot(3,1,1);
    plot(t, D.disp0, 'b-', 'LineWidth', 1.5); hold on;
    plot(t, D.disp_from_acc, 'r--', 'LineWidth', 1.5);
    xlabel('Time (s)', 'FontWeight', 'bold');
    ylabel('Displacement (m)', 'FontWeight', 'bold');
    legend('ANSYS disp (detrended)','Disp from acc (FFT integration)');
    title(sprintf('Broadband Displacement Reconstruction (corr=%.4f, NRMSE=%.2f%%)', ...
        M.broadband_corr, M.broadband_nrmse), 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    
    subplot(3,1,2);
    plot(t, C.disp_cardiac, 'b-', 'LineWidth', 1.5); hold on;
    plot(t, C.dispA_cardiac, 'r--', 'LineWidth', 1.5);
    xlabel('Time (s)', 'FontWeight', 'bold');
    ylabel('Displacement (m)', 'FontWeight', 'bold');
    legend('Disp cardiac (5-40 Hz)','Reconstructed disp cardiac');
    title(sprintf('Cardiac Band (5-40 Hz) - Waveform Match (corr=%.4f, NRMSE=%.2f%%)', ...
        M.cardiac_corr, M.cardiac_nrmse), 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    
    subplot(3,1,3);
    plot(t, R.disp_resp, 'b-', 'LineWidth', 1.5); hold on;
    plot(t, R.dispA_resp, 'r--', 'LineWidth', 1.5);
    plot(t, R.acc_env_resp, 'g:', 'LineWidth', 1.5);
    xlabel('Time (s)', 'FontWeight', 'bold');
    ylabel('Resp component (a.u.)', 'FontWeight', 'bold');
    legend('Resp from disp','Resp from disp(acc)','Resp from SCG envelope');
    title(sprintf('Respiration Band (0.1-0.5 Hz) - Agreement (corr=%.4f, lag=%.4fs)', ...
        M.resp_corr, M.resp_lagSec), 'FontSize', 14, 'FontWeight', 'bold');
    grid on;

    % Figure 2: Integration methods comparison
    figure('Name', 'Integration Methods Comparison', 'Position', [150, 150, 1400, 600]);
    plot(t, D.disp0, 'k-', 'LineWidth', 2); hold on;
    plot(t, D.disp_from_acc_polyfit, 'g--', 'LineWidth', 1.5);
    plot(t, D.disp_from_acc_fft, 'r-.', 'LineWidth', 1.5);
    xlabel('Time (s)', 'FontWeight', 'bold', 'FontSize', 12);
    ylabel('Displacement (m)', 'FontWeight', 'bold', 'FontSize', 12);
    legend('ANSYS disp','Polyfit detrend','FFT integration (used)', 'Location', 'best');
    title('Integration Method Comparison', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    
    % Figure 3: Sanity check - Acceleration comparison
    figure('Name', 'Sanity Check - Acceleration Consistency', 'Position', [200, 200, 1200, 600]);
    plot(t, D.acc0, 'b-', 'LineWidth', 1.5); hold on;
    plot(t, D.acc_from_disp, 'm--', 'LineWidth', 1.5);
    xlabel('Time (s)', 'FontWeight', 'bold', 'FontSize', 12);
    ylabel('Acceleration (m/s²)', 'FontWeight', 'bold', 'FontSize', 12);
    legend('ANSYS acc (actual)','Acc from disp (2x diff)', 'Location', 'best');
    title(sprintf('Sanity Check: Are signals from same node/axis? (corr=%.4f)', ...
        D.sanity_check_corr), 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    
    % Figure 4: Respiration rate tracking
    figure('Name', 'Respiration Rate Tracking', 'Position', [250, 250, 1200, 600]);
    plot(R.disp_rr.tmid, R.disp_rr.bpm, '-o', 'LineWidth', 2, 'MarkerSize', 8); hold on;
    plot(R.dispA_rr.tmid, R.dispA_rr.bpm, '-s', 'LineWidth', 2, 'MarkerSize', 8);
    plot(R.acc_env_rr.tmid, R.acc_env_rr.bpm, '-^', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('Time (s)', 'FontWeight', 'bold', 'FontSize', 12);
    ylabel('Resp rate (breaths/min)', 'FontWeight', 'bold', 'FontSize', 12);
    legend('RR from disp','RR from disp(acc)','RR from SCG envelope', 'Location', 'best');
    title(sprintf('Respiration Rate Tracking (MAE=%.2f bpm, RMSE=%.2f bpm)', ...
        M.rr_mae_bpm, M.rr_rmse_bpm), 'FontSize', 14, 'FontWeight', 'bold');
    grid on;

    % Figure 5: PSD snapshot example (first window)
    figure('Name', 'Power Spectral Density', 'Position', [300, 300, 1400, 800]);
    
    win = round(CFG.resp.winSec*fs);
    if win > length(R.disp_resp)
        win = length(R.disp_resp);
    end
    
    % Adjust Welch window size for short signals
    welch_win = min(round(win*0.5), win-1);
    welch_win = max(welch_win, 8);
    
    % PSD for displacement resp
    subplot(2,2,1);
    xw = R.disp_resp(1:win);
    [P,F] = pwelch(xw, hamming(welch_win), [], [], fs);
    plot(F,P, 'b-', 'LineWidth', 1.5);
    xlim([0 2]);
    xlabel('Frequency (Hz)', 'FontWeight', 'bold');
    ylabel('PSD', 'FontWeight', 'bold');
    title('PSD - Resp from Displacement', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    
    % PSD for reconstructed displacement resp
    subplot(2,2,2);
    xw = R.dispA_resp(1:win);
    [P,F] = pwelch(xw, hamming(welch_win), [], [], fs);
    plot(F,P, 'r-', 'LineWidth', 1.5);
    xlim([0 2]);
    xlabel('Frequency (Hz)', 'FontWeight', 'bold');
    ylabel('PSD', 'FontWeight', 'bold');
    title('PSD - Resp from Reconstructed Disp', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    
    % PSD for SCG envelope resp
    subplot(2,2,3);
    xw = R.acc_env_resp(1:win);
    [P,F] = pwelch(xw, hamming(welch_win), [], [], fs);
    plot(F,P, 'g-', 'LineWidth', 1.5);
    xlim([0 2]);
    xlabel('Frequency (Hz)', 'FontWeight', 'bold');
    ylabel('PSD', 'FontWeight', 'bold');
    title('PSD - Resp from SCG Envelope', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    
    % Overlay comparison
    subplot(2,2,4);
    xw1 = R.disp_resp(1:win);
    xw2 = R.dispA_resp(1:win);
    [P1,F1] = pwelch(xw1, hamming(welch_win), [], [], fs);
    [P2,F2] = pwelch(xw2, hamming(welch_win), [], [], fs);
    plot(F1,P1, 'b-', 'LineWidth', 1.5); hold on;
    plot(F2,P2, 'r--', 'LineWidth', 1.5);
    xlim([0 2]);
    xlabel('Frequency (Hz)', 'FontWeight', 'bold');
    ylabel('PSD', 'FontWeight', 'bold');
    title('PSD Comparison', 'FontSize', 12, 'FontWeight', 'bold');
    legend('Disp','Reconstructed Disp');
    grid on;
    
    % Figure 6: Cross-correlation
    figure('Name', 'Cross-Correlation Analysis', 'Position', [350, 350, 1000, 600]);
    maxLag = round(10*fs);
    [c,lags] = xcorr(zscore(R.disp_resp), zscore(R.dispA_resp), maxLag, 'coeff');
    plot(lags/fs, c, 'b-', 'LineWidth', 2);
    xlabel('Lag (seconds)', 'FontWeight', 'bold', 'FontSize', 12);
    ylabel('Correlation Coefficient', 'FontWeight', 'bold', 'FontSize', 12);
    title(sprintf('Cross-Correlation: Disp Resp vs Reconstructed Disp Resp\n(Max corr=%.4f at lag=%.4fs)', ...
        M.corr_disp_vs_dispA, M.lagSec), 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    hold on;
    plot(M.lagSec, M.corr_disp_vs_dispA, 'ro', 'MarkerSize', 12, 'LineWidth', 2);
    legend('Cross-correlation','Peak', 'Location', 'best');
end
