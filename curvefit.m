

%% New script for fitting a model to the raw data

% Model: y^ = A sin (wt+phi_0) + B

accelScale = 1/9.82; % scale accelerometer readings


M = readmatrix("recordings/recording_20250701_02.csv");
Gx = M(:,4);
Gy = M(:,5);
Gz = M(:,6);
G_mag = sqrt(Gx.^2+Gy.^2+Gz.^2);

epsilon = 1e-6;
correction_factor = ones(size(Gz));
valid = abs(Gz) > epsilon;

correction_factor(valid) = G_mag(valid) ./ abs(Gz(valid));
Gz_corrected  = correction_factor.*Gz;

% Lowpass 
fc = 10; 
fs = 100;
n = 100; % filter order
b = fir1(n, (fc/(fs/2)), 'low');

M_filt = M;
M_filt(:,1:3) = filtfilt(b,1,M(:,1:3));

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

    %if i == 1
    %    ylabel("Accel [g]")
    %elseif i == 4
    %    ylabel("Gyro [°/s]")
    %end

    if accelScale < 1
        ylabels = ["X[m/s²]", "Y [m/s²]", "Z [m/s²]", "X [°/s]", "Y [°/s]", "Z [°/s]"];
    else
        ylabels = ["X [g]", "Y [g]", "Z [g]", "X [°/s]", "Y [°/s]", "Z [°/s]"];
    end
    ylabel(ylabels(i))

    if i == 3 || i == 6
        xlabel("Time [s]")
    end
    
    %xtickformat("mm:ss")
end
linkaxes(ax,"x")

Ax = M_filt(:,1);
Ay = M_filt(:,2);
 

wheel_circ = 1.82;
angle_est = cumtrapz(-Gz)/100;

fprintf('Gyroscope distance estimate: %.2f m\n', angle_est(end)*wheel_circ/360);

%% DFT

fs = 100;
wsize = 100;
ovlap = wsize-1;
Ndft = 1024;

win = rectwin(wsize); % window function can be changed to something else
[sx,fx,tx, Px] = spectrogram(Ax,win,ovlap,Ndft,fs);
[sy,fy,ty, Py] = spectrogram(Ay,win,ovlap,Ndft,fs);


% car speed
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
gyro_vals = -Gz(wsize/2:end-wsize/2)*wheel_circ/100; % (DPS/360)*circ*3.6
p1 = plot(ty,gyro_vals,'Color',[1.0, 0.4, 0.0]);
legend(p1, 'Gyroscope signal overlay', 'Location', 'northwest');



%% Parameter estimation using DFT 

% b(0) = A
% b(1) = w
% b(2) = phi0
% b(3) = B
% yhat = A + sin(wt+phi0) + B

%model = @(t,b) b(1)*sin(2*pi*b(2)*t + b(3)) + b(4);

% b(2) is the total phase estimate theta_est.
model_dft = @(t,b) b{1}.*sin(b{2}) + b{3};

dt = 1/100;
N = size(M,1);   

% frequency and amplitude estimates
th = 1;   % threshold frequency in herz
th_idx = round(th/(fs/Ndft)) +1;   % index of roughly this frequency
threshold = (th_idx-1)*fs/Ndft;    % actual threshold with this index

[~,f0_relative_idx] = max(abs(sy(th_idx:end,:)));  % returns relative index of masked array
f0_idx = f0_relative_idx + th_idx - 1;            % corresponding index in full array

%non_dc_bins = 5:N/2;
%[~,peak_bin] = max(abs(sy(non_dc_bins)));
%freq_bin = non_dc_bins(peak_bin);

f_vals = fy(f0_idx);

% Prints speed estimate (not really needed)
f_avg = mean(fx(f0_idx(1600:5500)))
v_avg = f_avg*wheel_circ

cols = 1:size(sy,2);   % cols for sub2ind
Sx = sx(sub2ind(size(sx), f0_idx, cols));  % values of S at wanted frequencies
Sy = sy(sub2ind(size(sy), f0_idx, cols));


%A_est = (abs(Sy)+abs(Sx))/sum(win);
A_est = 2*abs(Sy/sum(win));
%A_est = 1;
% phase offset estimate (phi_0)
S = Sx + 1j*Sy;

tol = 1e-6;    % remove tiny noise
S(abs(S) < tol) = 0;

%theta_raw = unwrap(angle(S));
theta_raw = angle(S);
phi_offs = transpose(2*pi*f_vals*wsize/(2*fs));
theta_est = theta_raw+phi_offs;

% DC offset estimate
B_est = sy(1,:)/(sum(win));

b = {A_est, theta_est, B_est};
%b = {1, theta_est, B_est};

yhat = model_dft(ty,b);

figure;
dt = 1/100;
N = size(M,1);
t = (0:N-1)*dt;
t= transpose(t);

plot(t, M_filt(:,2)/accelScale)
hold on
plot(ty,yhat/accelScale)
grid on
%xtickformat("mm:ss")
xlabel("Time [s]");
if accelScale < 1
    ylabel("Acceleration [m/s²]");
else
    ylabel("Acceleration [g]");
end
% Then we could plot the unwrapped phase:
theta_unwr = unwrap(theta_raw);
theta_unwr = theta_unwr-theta_unwr(1); % remove initial phase
figure;
plot(ty, theta_unwr)
xlabel("Time [s]");
ylabel("Angle [rad]");

fprintf('STFT distance estimate: %.2f m\n', theta_unwr(end)*wheel_circ/(2*pi));

%% Parameter estimation using lsq 

window_size = 10;
N = length(Ay);
%step = 35;

model2 = @(b,t) b(1)*sin(2*pi*b(2)*t+b(3)) + b(4);

b0 = [0.1, 0.1, 0.0, 0]; % initial guess

lb = [0, 0, 0, min(Ay)];
ub = [4, 20, 2*pi, max(Ay)];

opts = optimoptions('lsqcurvefit', ...
    'MaxFunctionEvaluations', 5000, ...
     'Display', 'off'); % 'iter'

time = [];
vals = [];
params = {};

k=1;
for i=window_size/2:N-window_size/2
    j = i-window_size/2+1:(i+window_size/2);
    t_win = t(j);
    Ay_win = Ay(j);
    
    t_rel = t_win-t_win(1);

    b_hat = lsqcurvefit(model2,b0,t_rel,Ay_win, lb, ub, opts);
    y_fit = model2(b_hat,t_rel);

    time(k) = t_win(window_size/2);
    vals(k) = y_fit(window_size/2);
    params{k} = b_hat;

    b0 = b_hat;
    k=k+1;
end

figure;
%colors = lines(10); 
%num_colors = size(colors, 1);


plot(t,Ay/accelScale, 'Color','k')
hold on
    %c = colors(mod(i-1, num_colors) + 1, :);  
plot(time,vals/accelScale,'Color', 'r')
if accelScale < 1
    ylabel('Acceleration [m/s²]')
else
    ylabel('Acceleration [g]')
end
xlabel('Time [s]')

%% Distance calculation

param_matrix = cell2mat(params');
freqs = param_matrix(:,2);
phis = param_matrix(:,3);

%t_offset = (window_size/2) / fs;
%t_aligned = time(:)- t_offset;
%t_stripped = t_aligned(3:end-2);
%time = time(:);

phase = 2*pi*freqs.*time(:)+phis;
phase_unwr = unwrap(phase);

figure;
plot(phase_unwr)

fprintf('LSQ distance estimate: %.2f\n', theta_unwr(end)*wheel_circ/(2*pi));
