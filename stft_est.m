

%%  STFT estimation

clear

accelScale = 1/9.82;

M = readmatrix("recordings/recording_20250701_04.csv");
Gx = M(:,4);
Gy = M(:,5);
Gz = M(:,6);

% Lowpass 
fc = 6; 
fs = 100;
n = 100; % filter order
b = fir1(n, (fc/(fs/2)), 'low');

M_filt = M;
M_filt(:,1:3) = filtfilt(b,1,M(:,1:3));

Ax = M_filt(:,1);
Ay = M_filt(:,2);

%%
figure("Name","raw+filtered")
dt = 1/100;
N = size(M,1);
t = (0:N-1)*dt;
t= transpose(t);
tiledlayout(3,2,"TileIndexing","columnmajor")
ax = [];
for i = 1:6
    ax(end+1) =nexttile;
    if i<4
        plot(t, M(:,i)/accelScale)
        hold on
        plot(t,M_filt(:,i)/accelScale)
    else
        plot(t, M(:,i))
    end

    if i == 6
        plot(t,Gz)
    end
    grid on

    if accelScale < 1
        ylabels = ["X[m/s²]", "Y [m/s²]", "Z [m/s²]", "X [°/s]", "Y [°/s]", "Z [°/s]"];
    else
        ylabels = ["X [g]", "Y [g]", "Z [g]", "X [°/s]", "Y [°/s]", "Z [°/s]"];
    end
    ylabel(ylabels(i))

    if i == 3 || i == 6
        xlabel("Time [s]")
    end
end
linkaxes(ax,"x")


%% STFT stuff

fs = 100;
wsize = 130;
ovlap = wsize-1;
Ndft = 1024;

%win = hann(wsize); % window function can be changed to something else
win = rectwin(wsize);
[sx,fx,tx, Px] = spectrogram(Ax,win,ovlap,Ndft,fs);
[sy,fy,ty, Py] = spectrogram(Ay,win,ovlap,Ndft,fs);

%%
wheel_circ = 1.82;
v_car = fy*wheel_circ*3.6;

figure;
imagesc(ty, v_car, 10*log10(Py));
axis xy;
xlabel("Time [s]");
ylabel("Speed [km/h]");
title("Spectrogram with overlay");
grid on;
c = colorbar;
c.Label.String = 'Power/frequency [dB/Hz]';

hold on;
gyro_vals = -Gz(wsize/2:end-wsize/2)*wheel_circ/100; % (DPS/360)*circ*3.6 [km/h]
p1 = plot(ty,gyro_vals,'Color',[1.0, 0.4, 0.0]);
legend(p1, 'Gyroscope signal overlay', 'Location', 'northwest');


%% Report plots

win_func = @hann;


start_indices = [270, 290,310,330, 340,350];
%start_indices = 290;
segment = Ax(start_indices:start_indices + wsize - 1);

rect_win = ones(wsize,1);
rect_segment = segment .* rect_win;

win = win_func(wsize);
tapered_segment = segment .* win;

Nfft = 1024;
f = (0:Nfft-1) * (fs/Nfft);
t = (0:wsize-1)/fs;

X_rect = fft(rect_segment, Nfft);
X_tapered = fft(tapered_segment, Nfft);

X_rect_mag = abs(X_rect);
X_tapered_mag = abs(X_tapered);


%%
figure;
subplot(1,2,1);
t = (0:wsize-1)/fs;
plot(t, rect_segment, 'k-', 'DisplayName', 'Raw segment');
hold on;
plot(t, tapered_segment, 'r-', 'DisplayName', 'Windowed');
xlabel('Time [s]');
ylabel('Amplitude');
title('Time-domain Window');
legend;
grid on;

% Frequency-domain (FFT)
subplot(1,2,2);
%plot(f, X_rect_mag, 'b', 'DisplayName','rect');
hold on
plot(f,X_tapered_mag);
xlim([0 fs/2]);
xlabel('Frequency [Hz]');
ylabel('|X(f)|');
title('FFT of Windowed Signal');
%legend;
grid on;


%%


%start_indices = [270, 280, 290];
fig = figure('Units','normalized','OuterPosition',[0 0 1 1]); 
set(fig, 'PaperOrientation', 'landscape');

%figure;
subplot(1,2,1);
hold on;
xlabel('Time [s]');
ylabel('Amplitude');
title('Time-domain Windows');
grid on;

subplot(1,2,2);
hold on;
xlim([0 fs/2]);
xlabel('Frequency [Hz]');
ylabel('|X(f)|');
title('FFT of Windowed Segments');
grid on;


for i = 1:length(start_indices)
    idx = start_indices(i);
    segment = Ax(idx : idx + wsize - 1);

    % Apply windows
    win = win_func(wsize);
    tapered_segment = segment .* win;

    % FFT
    X_tapered = fft(tapered_segment, Nfft);
    X_mag = abs(X_tapered);

    % Plot time-domain signal
    subplot(1,2,1);
    plot(t + idx/fs, tapered_segment, 'DisplayName', sprintf('Frame %d', i)); % offset time axis

    % Plot frequency-domain magnitude
    subplot(1,2,2);
    plot(f, X_mag, 'DisplayName', sprintf('Frame %d', i));
end

subplot(1,2,1);
t2 = (start_indices(1):(start_indices(end)+wsize-1))/fs;
plot(t2, Ax(start_indices(1):start_indices(end) + wsize - 1), 'DisplayName','Raw segment',LineWidth=1.5)


% Add legends after loop
subplot(1,2,1); legend;
subplot(1,2,2); legend;
%%
%print(fig, 'my_figure.pdf', '-dpdf', '-fillpage');
fig = gcf;
exportgraphics(fig, 'FFTFrames_01.pdf', 'ContentType', 'vector');


%% Parameter estimation using STFT 

model_dft = @(t,b) cos(b{1}) + b{2};

dt = 1/100;
N = size(M,1);   

% 1 Hz = 1.82 m/s = 6.55 km/h
%th = 1.5;   % threshold frequency in Hz
%th_idx = round(th/(fs/Ndft)) +1;   % index of roughly this frequency
%threshold = (th_idx-1)*fs/Ndft;    % actual threshold with this index

%[~,f0y_relative_idx] = max(abs(sy(th_idx:end,:)));  % returns relative index of masked array
%f0_idx = f0y_relative_idx + th_idx - 1;            % corresponding index in full array

%DC_amplitude = abs(sy(1,:));
%target_amplitude = max(abs(sy(th_idx:end,:)),[],1);
%motion_detected = target_amplitude > 21.5; % This value determines when we start looking only above the threshold.

% NEW threshold logic
gyro_vals = -Gz(wsize/2:end-wsize/2)/360; % DPS/360 [1/s]
freq_res = fs/Ndft;

f_th = 1.7; % Roughly 6 km/h
th_idx = round(f_th / freq_res) + 1;

f0_idx = nan(1,size(sy,2));
%%
for t = 1:size(sy,2)
    if gyro_vals(t)<f_th
        % use gyro derived freq
        idx = round(gyro_vals(t) / freq_res) + 1;
        idx = max(1, min(idx, Ndft));
        f0_idx(t) = idx;
    else
        % use STFT peak
        [~,rel_idx] = max(abs(sy(th_idx:end,t)));
        idx = rel_idx + th_idx - 1;
        f0_idx(t) = idx;
    end
end

%%

%dc_magnitude = abs(sy(1,:));
%dc_zero_idx = find(dc_magnitude<1, 1 , 'first'); % index in time where dc is 0

% Up until dc_zero_idx
%f0_idx_dyn = zeros(1,size(sy,2));
%[~, f0_idx_dyn(1:dc_zero_idx)] = max(abs(sy(:, 1:dc_zero_idx)), [], 1);


% [~, temp_idx] = max(abs(sy(15:end, dc_zero_idx+1:end)), [], 1);
% f0_idx_dyn(dc_zero_idx+1:end) = temp_idx + 14;

% Then bandpass
Py_bp = zeros(size(sy));
bp_width = 0.5; % in Hz
curr_max_idx = f0_idx_dyn(dc_zero_idx); % current maximum frequency peak index
for t_idx = dc_zero_idx+1:size(sy,2)
        
    % obtain bp range from prev. window
    lo_f = max(0, fy(curr_max_idx) - bp_width);
    hi_f = fy(curr_max_idx) + bp_width;

    f_pass = fy >= lo_f & fy <= hi_f;
    sy_bp = zeros(size(sy(:,t_idx)));
    sy_bp(f_pass) = sy(f_pass, t_idx);
    
    Py_bp(:, t_idx) = abs(sy_bp).^2 / (norm(win)^2);

    [~,max_idx] = max(abs(sy_bp));
    f0_idx_dyn(t_idx) = max_idx; % store
    curr_max_idx = max_idx; % update
end

%%

figure;
plot(ty,dc_magnitude);
hold on
dc_magnitude = abs(sy(1,:)) + abs(sy(2,:)) + abs(sy(3,:));
plot(ty,dc_magnitude);

%%
plot(ty, f0_idx, 'DisplayName','f0_idx')
hold on 
plot(ty, f0_idx_dyn, 'DisplayName','f0_idx_dyn')
legend
%%

f_peak = fy(f0_idx_dyn);        % Hz
v_peak = f_peak * wheel_circ * 3.6;

v_lo = max(v_peak - bp_width*wheel_circ*3.6, 0);
v_hi = v_peak + bp_width*wheel_circ*3.6;

figure;
imagesc(ty, v_car, 10*log10(Py_bp));
axis xy;
xlabel('Time (s)');
ylabel('Speed (km/h)');
title('Spectrogram');
colorbar;

hold on;
plot(ty, v_lo, 'w--', 'LineWidth', 1.5);  % Lower band edge
plot(ty, v_hi, 'w--', 'LineWidth', 1.5);  % Upper band edge
plot(ty, v_peak, 'r-', 'LineWidth', 2);   % Tracked speed
legend('Lower Band Edge', 'Upper Band Edge', 'Tracked Speed');


%%

% f0_idx = nan(1,size(sy,2));
% for t = 1:size(sy,2)
%     if motion_detected(t)
%         [~,relative_idx] = max(abs(sy(th_idx:end,(t))));
%         f0_idx(t) = relative_idx + th_idx - 1; 
%     else
%         f0_idx(t) = 1;
%     end
% end
%%
f_vals = fy(f0_idx_dyn);

cols = 1:size(sy,2);   % cols for sub2ind
Sx = sx(sub2ind(size(sx), f0_idx_dyn, cols));  % values of S at wanted frequencies
Sy = sy(sub2ind(size(sy), f0_idx_dyn, cols));

% phase offset estimate (phi_0)
S = Sx + 1j*Sy;

tol = 1e-6;    % remove tiny noise
S(abs(S) < tol) = 0;

%theta_dft = unwrap(angle(S));
theta_dft = angle(S);
phiy_offs = transpose(2*pi*f_vals*wsize/(2*fs));
thetay_est = theta_dft+phiy_offs;

phix_offs = transpose(2*pi*f_vals*wsize/(2*fs));
thetax_est = theta_dft+phix_offs;

% DC offset estimate
% By_est = sy(1,:)/(sum(win));
% Bx_est = sx(1,:)/(sum(win));
% dc_frames = (gyro_vals<f_th);
% alphax = -0.01;
% alphay = 0.005;
% By_est(dc_frames) = 2*pi*alphay*gyro_vals(dc_frames).^2;
% Bx_est(dc_frames) = 2*pi*alphax*gyro_vals(dc_frames).^2;

by = {thetay_est-pi/2, 0};
bx = {thetax_est, 0};

yhat = model_dft(ty,by);
xhat = model_dft(tx,bx);

%yhat_corrected = yhat-By_est;
%xhat_corrected = xhat-Bx_est;

figure('Name','Raw + reconstructed');
dt = 1/100;
N = size(M,1);
t = (0:N-1)*dt;
t= transpose(t);
tl = tiledlayout(2,1,"TileIndexing","columnmajor");
xlabel(tl,'Time [s]');

if accelScale < 1
    ylabel(tl,"Acceleration [m/s²]");
else
    ylabel(tl,"Acceleration [g]");
end

ax = [];
ax(end+1) =nexttile;
plot(t, M_filt(:,1)/accelScale,'DisplayName','raw x')
hold on
plot(tx,xhat/accelScale,'DisplayName','xhat')

p = 20;
for i = 1:p:length(ty)
    text(tx(i), xhat(i), ...
         sprintf('%.2f', f_vals(i)), ...
         'FontSize', 6, ...
         'Color', 'blue', ...
         'VerticalAlignment', 'bottom', ...
         'HorizontalAlignment', 'center');
end

title('x-axis')
grid on
legend

ax(end+i)=nexttile;
plot(t, M_filt(:,2)/accelScale,'DisplayName','raw y')
hold on
plot(ty,yhat/accelScale,'DisplayName','yhat')

for i = 1:p:length(ty)
    text(ty(i), yhat(i), ...
         sprintf('%.2f', f_vals(i)), ...
         'FontSize', 6, ...
         'Color', 'blue', ...
         'VerticalAlignment', 'bottom', ...
         'HorizontalAlignment', 'center');
end
title('y-axis')
grid on
legend

linkaxes(ax,"x")

%% 
figure;
plot(xhat, yhat)
