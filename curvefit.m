

%% New script for fitting a model to the raw data

% Model: y^ = A sin (wt+phi_0) + B

accelScale = 1/9.82; % scale accelerometer readings


M = readmatrix("recordings/recording_20250701_03.csv");
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
fc = 6; 
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
model_dft = @(t,b) cos(b{1}) + b{2};

dt = 1/100;
N = size(M,1);   

% 1 Hz = 1.82 m/s = 6.55 km/h
% frequency and amplitude estimates
th = 1;   % threshold frequency in Hz
th_idx = round(th/(fs/Ndft)) +1;   % index of roughly this frequency
threshold = (th_idx-1)*fs/Ndft;    % actual threshold with this index

[~,f0y_relative_idx] = max(abs(sy(th_idx:end,:)));  % returns relative index of masked array
%[~,f0x_relative_idx] = max(abs(sx(th_idx:end,:)));  % returns relative index of masked array

f0_idx = f0y_relative_idx + th_idx - 1;            % corresponding index in full array
%f0x_idx = f0x_relative_idx + th_idx - 1;            % corresponding index in full array
%f0_idx = round((f0y_idx+f0x_idx)/2);

DC_amplitude = abs(sy(1,:));
target_amplitude = max(abs(sy(th_idx:end,:)),[],1);
%ratio = target_amplitude ./ DC_amplitude;


motion_detected = target_amplitude > 25; % This value determines when we start looking only above the threshold. 

f0_idx = nan(1,size(sy,2));

for t = 1:size(sy,2)
    if motion_detected(t)
        [~,relative_idx] = max(abs(sy(th_idx:end,(t))));
        f0_idx(t) = relative_idx + th_idx - 1; 
    else
        f0_idx(t) = 1;
    end
end

f_vals = fy(f0_idx);

% Prints speed estimate (not really needed)
f_avg = mean(fx(f0_idx(1600:5500)))
v_avg = f_avg*wheel_circ

cols = 1:size(sy,2);   % cols for sub2ind
Sx = sx(sub2ind(size(sx), f0_idx, cols));  % values of S at wanted frequencies
Sy = sy(sub2ind(size(sy), f0_idx, cols));

% phase offset estimate (phi_0)
S = Sx + 1j*Sy;

tol = 1e-6;    % remove tiny noise
S(abs(S) < tol) = 0;

%theta_raw = unwrap(angle(S));
theta_dft = angle(S);
phiy_offs = transpose(2*pi*f_vals*wsize/(2*fs));
thetay_est = theta_dft+phiy_offs;

phix_offs = transpose(2*pi*f_vals*wsize/(2*fs));
thetax_est = theta_dft+phix_offs;

% DC offset estimate
By_est = sy(1,:)/(sum(win));
Bx_est = sx(1,:)/(sum(win));
%Bx_est(f0_idx==1) = 0;
%By_est(f0_idx==1) = 0;
%Bx_est=0;
%By_est=0;

by = {thetay_est-pi/2, By_est};
bx = {thetax_est, Bx_est};

% Want this to work
%by = {1, thetay_est, By_est};
%bx = {1, thetax_est, Bx_est};

%b = {1, theta_est, B_est};

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
title('x-axis')
grid on
legend

ax(end+i)=nexttile;
plot(t, M_filt(:,2)/accelScale,'DisplayName','raw y')
hold on
plot(ty,yhat/accelScale,'DisplayName','yhat')
title('y-axis')
grid on
legend

linkaxes(ax,"x")

%%

% Then we could plot the unwrapped phase:
%theta_dft = theta_dft;
%theta_dft = unwrap(theta_dft)
theta_dft_wrapped = wrapToPi(thetay_est - thetay_est(1));

theta_alt = atan2(yhat_corrected,xhat_corrected);
theta_alt = wrapToPi(theta_alt - theta_alt(1));
%theta_alt = unwrap(theta_alt-theta_alt(1));

theta_raw= atan2(Ay,Ax);
theta_raw = wrapToPi(theta_raw - theta_raw(1));

%theta_raw = theta_raw;
%theta_raw = unwrap(theta_raw-theta_raw(1));

figure;
plot(ty, theta_dft_wrapped, 'DisplayName','Angle(S)')
hold on
plot(ty,theta_alt+0.01, 'DisplayName', 'atan2(yhat,xhat)')
plot(t,theta_raw,'DisplayName','arg raw')
xlabel("Time [s]");
ylabel("Angle [rad]");
legend
grid on

%%

Ay_corrected = Ay(wsize/2:end-wsize/2) - By_est.';
Ax_corrected = Ax(wsize/2:end-wsize/2) - Bx_est.';

figure;
plot(Ax,Ay, 'DisplayName','raw')
hold on
plot(Ax_corrected, Ay_corrected, 'DisplayName', 'corrected')
xlabel('x')
ylabel('y')
legend

fprintf('STFT distance estimate: %.2f m\n', theta_dft(end)*wheel_circ/(2*pi));

%% Parameter estimation using lsq 

window_size = 10;
N = length(Ay);

model = @(b,t) sin(2*pi*b(1)*t+b(2)) + b(3);
b0 = [0.1, 0.0, 0.0]; % initial guess

lb = [-5.5, -pi, min(Ay)];
ub = [5.5, pi, max(Ay)];

opts = optimoptions('lsqcurvefit', ...
    'MaxFunctionEvaluations', 5000, ...
     'Display', 'off'); % 'iter'

time = [];
xvals = [];
yvals = [];
xparams = {};
yparams = {};

k=1;
for i=window_size/2:4:N-window_size/2
    j = i-window_size/2+1:(i+window_size/2);
    t_win = t(j);
    Ay_win = Ay(j);
    Ax_win = Ax(j);
    
    t_rel = t_win-t_win(1);

    b_hatx = lsqcurvefit(model, b0, t_rel, Ax_win, lb, ub, opts);
    b_haty = lsqcurvefit(model, b0, t_rel, Ay_win, lb, ub, opts);

    y_fit = model(b_haty,t_rel);
    x_fit = model(b_hatx,t_rel);

    time(k) = t_win(window_size/2);
    yvals(k) = y_fit(window_size/2);
    xvals(k) = x_fit(window_size/2);
    xparams{k} = b_hatx;
    yparams{k} = b_haty;

    %b0 = b_hat;
    k=k+1;
end


%%
figure;

plot(t,Ay/accelScale, 'Color','k')
hold on
plot(time,yvals/accelScale,'Color', 'r')

if accelScale < 1
    ylabel('Acceleration [m/s²]')
else
    ylabel('Acceleration [g]')
end
xlabel('Time [s]')
title('Y-axis')

figure;

plot(t,Ax/accelScale, 'Color','k')
hold on
plot(time,xvals/accelScale,'Color', 'r')

if accelScale < 1
    ylabel('Acceleration [m/s²]')
else
    ylabel('Acceleration [g]')
end
xlabel('Time [s]')
title('X-axis')

%% Distance calculation

xparam_matrix = cell2mat(xparams');
yparam_matrix = cell2mat(yparams');

xoffs = xparam_matrix(:,3);
yoffs = yparam_matrix(:,3);

xcorr = xvals - xoffs.';
ycorr = yvals - yoffs.';

figure;
plot(time,xvals)
hold on
plot(time,xcorr)
plot(time,yvals)

figure;
plot(xvals,yvals)
hold on
plot(xcorr,ycorr)
%%

arg_x = xparam_matrix(:,2);
arg_y = yparam_matrix(:,2);

phase_raw = atan2(yvals,xvals);
phase_corr = atan2(ycorr,xcorr);
phase_unwr = unwrap(phase_raw);

figure;
plot(time,phase_corr,'DisplayName','phase (offset removed)')
hold on
plot(time,phase_raw,'DisplayName','raw phase')
plot(time,arg_x,'DisplayName','omega*t+phi')
ylabel('Angle [rad]')
xlabel('Time [s]')
title('LSQ Angle Estimate')
legend

fprintf('LSQ distance estimate: %.2f\n', phase_unwr(end)*wheel_circ/(2*pi));