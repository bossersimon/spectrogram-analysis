

%%  STFT estimation

accelScale = 1;%/9.82;

M = readmatrix("recordings/recording_20250701_02.csv");
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
        plot(t,Gz_corrected)
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
wsize = 100;
ovlap = wsize-1;
Ndft = 4096;

win = hann(wsize); % window function can be changed to something else
[sx,fx,tx, Px] = spectrogram(Ax,win,ovlap,Ndft,fs);
[sy,fy,ty, Py] = spectrogram(Ay,win,ovlap,Ndft,fs);

%%
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


f_vals = fy(f0_idx);

cols = 1:size(sy,2);   % cols for sub2ind
Sx = sx(sub2ind(size(sx), f0_idx, cols));  % values of S at wanted frequencies
Sy = sy(sub2ind(size(sy), f0_idx, cols));

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
By_est = sy(1,:)/(sum(win));
Bx_est = sx(1,:)/(sum(win));
dc_frames = (gyro_vals<f_th);
alphax = -0.01;
alphay = 0.005;
By_est(dc_frames) = 2*pi*alphay*gyro_vals(dc_frames).^2;
Bx_est(dc_frames) = 2*pi*alphax*gyro_vals(dc_frames).^2;

by = {thetay_est-pi/2, By_est};
bx = {thetax_est, Bx_est};

yhat = model_dft(ty,by);
xhat = model_dft(tx,bx);

yhat_corrected = yhat-By_est;
xhat_corrected = xhat-Bx_est;

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
