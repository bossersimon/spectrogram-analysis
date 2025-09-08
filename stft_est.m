
%%  STFT estimation

clear

accelScale = 1/9.82;

M = readmatrix("recordings/recording_20250701_06.csv");
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
win = gausswin(wsize);
%win = rectwin(wsize);
[sx,fx,tx, Px] = spectrogram(Ax,win,ovlap,Ndft,fs);
[sy,fy,ty, Py] = spectrogram(Ay,win,ovlap,Ndft,fs);

wheel_circ = 1.82;
v_car = fy*wheel_circ*3.6; % speed in km/h

%%

fig = figure('Units','normalized','OuterPosition',[0 0 1 1]); 
set(fig, 'PaperOrientation', 'landscape');

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
legend(p1, 'Gyroscope signal overlay', 'Location', 'northwest','FontSize', 18);


%% Parameter estimation using STFT 

dt = 1/100;
N = size(M,1);   

% NEW threshold logic
gyro_vals = -Gz(wsize/2:end-wsize/2)/360; % DPS/360 [1/s]
freq_res = fs/Ndft;

f_th = 1; % Roughly 6 km/h
th_idx = round(f_th / freq_res) + 1;

f0_idx = nan(1,size(sy,2));

%%

curr_max_idx = 1;

%Py_bp = zeros(size(sy));
bp_width = 0.3; % in Hz
max_jump = 3;
for t = 1:size(sy,2)
    if gyro_vals(t)<f_th 
        % use gyro derived freq
        idx = round(gyro_vals(t) / freq_res) + 1;
        idx = max(1, min(idx, Ndft));

        if abs(idx - curr_max_idx) <= max_jump
            f0_idx(t) = idx;
            curr_max_idx = idx;
        else
            f0_idx(t) = curr_max_idx;  % ignore jump
        end
    else
        % obtain bp range from prev. window
        lo_f = max(0, fy(curr_max_idx) - bp_width);
        hi_f = fy(curr_max_idx) + bp_width;
    
        f_pass = fy >= lo_f & fy <= hi_f;
        sy_bp = zeros(size(sy(:,t)));
        sy_bp(f_pass) = sy(f_pass, t);
        
        %Py_bp(:, t) = abs(sy_bp).^2 / (norm(win)^2);
    
        [~,max_idx] = max(abs(sy_bp));
        f0_idx(t) = max_idx; % store
        curr_max_idx = max_idx; % update
    end
end


%%

f0_idx_dyn = zeros(1,size(sy,2));

Py_bp = zeros(size(sy));
bp_width = 0.3; % in Hz
curr_max_idx = 1; % current maximum frequency peak index
for t_idx = 1:size(sy,2)
        
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
dc_magnitude = abs(sy(1,:));
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

f_peak = fy(f0_idx);        % Hz
v_peak = f_peak * wheel_circ * 3.6;

v_lo = max(v_peak - bp_width*wheel_circ*3.6, 0);
v_hi = v_peak + bp_width*wheel_circ*3.6;

f_peak2 = fy(f0_idx_dyn);        % Hz
v_peak2 = f_peak2 * wheel_circ * 3.6;

v_lo2 = max(v_peak2 - bp_width*wheel_circ*3.6, 0);
v_hi2 = v_peak2 + bp_width*wheel_circ*3.6;

fig = figure('Units','normalized','OuterPosition',[0 0 1 1]); 
set(fig, 'PaperOrientation', 'landscape');

imagesc(ty, v_car, 10*log10(Py));
axis xy;
xlabel('Time (s)');
ylabel('Speed (km/h)');
title('Spectrogram');
colormap(bone);
colorbar;

colors = {'k', [0.3 0.3 0.3], [0.6 0.6 0.6], [0.9 0.9 0.9]};
styles = {'-', '--', '-.', ':'};

hold on;
plot(ty, v_lo, 'Color', colors{1},'LineStyle', styles{1}, 'LineWidth', 1.5);  % Lower band edge
plot(ty, v_hi, 'b--', 'LineWidth', 1.5);  % Upper band edge
plot(ty, v_peak, 'r-', 'LineWidth', 2);   % Tracked speed

plot(ty, v_lo2, 'w--', 'LineWidth', 1.5);  % Lower band edge
plot(ty, v_hi2, 'w--', 'LineWidth', 1.5);  % Upper band edge
plot(ty, v_peak2, 'r-', 'LineWidth', 2);   % Tracked speed
legend('Lower Band Edge', 'Upper Band Edge', 'Tracked Speed','FontSize',18);


%%
f_vals = fy(f0_idx);

model_dft = @(t,b) cos(b{1}) + b{2};

cols = 1:size(sy,2);   % cols for sub2ind
Sx = sx(sub2ind(size(sx), f0_idx, cols));  % values of S at wanted frequencies
Sy = sy(sub2ind(size(sy), f0_idx, cols));

% phase offset estimate (phi_0)
S = Sx + 1j*Sy;

tol = 1e-6;    % remove tiny noise
S(abs(S) < tol) = 0;

%theta_dft = unwrap(angle(S));
theta_dft = angle(S);

phi_offs = transpose(2*pi*f_vals*wsize/(2*fs));
thetay_est = theta_dft+phi_offs;
thetax_est = theta_dft+phi_offs;

by = {thetay_est-pi/2, 0};
bx = {thetax_est, 0};

yhat = model_dft(ty,by);
xhat = model_dft(tx,bx);

theta2 = atan2(yhat,xhat);

ax = [];

figure('Name','Raw + reconstructed');
dt = 1/100;
N = size(M,1);
t = (0:N-1)*dt;
t= transpose(t);
tl = tiledlayout(2,1,"TileIndexing","columnmajor");
xlabel(tl,'Time [s]','FontSize', 20, 'Interpreter', 'latex');

if accelScale < 1
    ylabel(tl,"Acceleration [m/s$^2$]", 'Interpreter','latex', 'FontSize', 18);
else
    ylabel(tl,"Acceleration [g]", 'Interpreter','latex', 'FontSize', 18);
end

%ax = [];
ax(end+1) =nexttile;
plot(t, Ax/accelScale,'DisplayName','raw x','LineWidth',1.5)
hold on
plot(tx,xhat/accelScale,'DisplayName','xhat','LineWidth',1.5)

title('x-axis','FontSize', 20, 'Interpreter', 'latex')
grid on
legend('FontSize', 16)

ax(end+i)=nexttile;
plot(t, Ay/accelScale,'DisplayName','raw y','LineWidth',1.5)
hold on
plot(ty,yhat/accelScale,'DisplayName','yhat','LineWidth',1.5)

title('y-axis','FontSize', 20, 'Interpreter', 'latex')
grid on
legend('FontSize', 16)

linkaxes(ax,"x")

%%

dphi = diff(thetax_est);

plot(ty(1:end-1), dphi)


%%
fig = gcf;
exportgraphics(fig, 'raw+recon_03.pdf', 'ContentType', 'vector');

%%


fig = figure('Units','normalized','OuterPosition',[0 0 1 1]); 
set(fig, 'PaperOrientation', 'landscape');

dt = 1/100;
N = size(M,1);
t = (0:N-1)*dt;
t= transpose(t);
tl = tiledlayout(2,1,"TileIndexing","columnmajor");
xlabel(tl,'Time [s]','Interpreter','latex', 'FontSize', 20);

ax = [];
ax(end+1) =nexttile;
plot(ty,yhat/accelScale,'DisplayName','yhat', 'LineWidth',1.5)

hold on
plot(tx,xhat/accelScale,'DisplayName','xhat','LineWidth',1.5)

if accelScale < 1
    ylabel("Acceleration [m/s$^2$]", 'Interpreter','latex', 'FontSize', 18);
else
    ylabel("Acceleration [g]", 'Interpreter','latex', 'FontSize', 18);
end

title('Time signals', 'FontSize', 20, 'Interpreter', 'latex');
grid on
legend

ax(end+i)=nexttile;
tht = wrapToPi(thetax_est);

plot(ty, tht,'DisplayName','phase','LineWidth',1.5,'Color','b')
hold on
ylabel("Phase [rad]", 'FontSize', 18, 'Interpreter','latex');
title('Phase', 'FontSize', 20, 'Interpreter', 'latex');
grid on

linkaxes(ax,"x")
%%
plot(ty,tht)


%% export

fig = gcf;
exportgraphics(fig, 'phasedist_02.pdf', 'ContentType', 'vector');


%% 
figure;
plot(xhat, yhat)


%%

fc = 1; 
fs = 100;
n = 100; % filter order
b = fir1(n, (fc/(fs/2)), 'high');


Ax_filtered = filtfilt(b,1,Ax(wsize/2:end-wsize/2));
Ay_filtered = filtfilt(b,1,Ay(wsize/2:end-wsize/2));
xhat_filtered = filtfilt(b,1,xhat);

figure;
plot(ty,Ax_filtered)
hold on
plot(ty,xhat_filtered)


mse = sum((Ax_filtered(:)-xhat(:)).^2)/(size(xhat,2));
%%
theta_raw  = atan2(Ay_filtered,Ax_filtered);
theta_raw2 = atan2(Ay,Ax);


fig = figure('Units','normalized','OuterPosition',[0 0 1 1]); 
set(fig, 'PaperOrientation', 'landscape');

%plot(ty,theta_raw, 'DisplayName','raw no offs')
hold on
%plot(ty,tht, 'DisplayName','est')
plot(t,theta_raw2, 'DisplayName','raw','Color','b', 'LineWidth',1.5)
title("raw phase")
xlabel('Time [s]')
ylabel('Angle [rad]')
%legend

%% new figure


fig = figure('Units','normalized','OuterPosition',[0 0 1 1]); 
set(fig, 'PaperOrientation', 'landscape');

dt = 1/100;
N = size(M,1);
t = (0:N-1)*dt;
t= transpose(t);
tl = tiledlayout(2,2,"TileIndexing","rowmajor");
xlabel(tl,'Time [s]','FontSize', 20, 'Interpreter', 'latex');

if accelScale < 1
    ylabel(tl,"Acceleration [m/s$^2$]", 'Interpreter','latex', 'FontSize', 18);
else
    ylabel(tl,"Acceleration [g]", 'Interpreter','latex', 'FontSize', 18);
end

ax = [];
ax(end+1) =nexttile(1);
plot(t, Ax/accelScale,'DisplayName','raw x','LineWidth',1.5)
hold on
plot(tx,xhat/accelScale,'DisplayName','xhat', 'LineWidth',1.5)

title('x-axis','FontSize', 20, 'Interpreter', 'latex')
grid on
legend('FontSize', 16)

ax(end+1)=nexttile(3);
plot(t, Ay/accelScale,'DisplayName','raw y','LineWidth',1.5)
hold on
plot(ty,yhat/accelScale,'DisplayName','yhat','LineWidth',1.5)

title('y-axis','FontSize', 20, 'Interpreter', 'latex')
grid on
legend('FontSize', 16)


ax(end+1) = nexttile(2, [2 1]);
plot(ty, tht, 'DisplayName','phi_{est}')
hold on
plot(t,theta_raw2, 'DisplayName', 'raw phase','LineWidth',1.5)

title('phase','FontSize', 20, 'Interpreter', 'latex','LineWidth',1.5)
grid on
legend('FontSize', 16)
ylabel("Angle [rad]", 'Interpreter','latex', 'FontSize', 18);


linkaxes(ax,"x")


%% Distance, speed calculation


gpslog = readmatrix("recordings/GPS/sensorLog_20250701_06.txt");

time= gpslog(:,1);
lat = gpslog(:,3);
lon = gpslog(:,4);
alt = gpslog(:,5);


fig = figure('Units','normalized','OuterPosition',[0 0 1 1]); 
set(fig, 'PaperOrientation', 'landscape');

geoplot(lat,lon,"-.",'LineWidth',1,'Color','r')
geobasemap streets
title("Car route")



lat0 = lat(1);
lon0 = lon(1);
alt0 = alt(1);

ref = [lat0 lon0 alt0];

wgs84 = wgs84Ellipsoid('meters');
[x, y, zUp] = geodetic2enu(lat, lon, alt, lat0, lon0, alt0, wgs84);

figure;
plot(x, y, 'b.-');
hold on
plot(x(1),y(1),'r*','DisplayName','Start')
plot(x(end),y(end),'g*', 'DisplayName', 'End')
axis equal;
xlabel('East [m]');
ylabel('North [m]');
legend

dx = diff(x);
dy = diff(y);

time_sec = time/1000;
datetime_array = datetime(time_sec, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');
sample_times = seconds(datetime_array - datetime_array(1));
disp(sample_times);

distances = (sqrt(dx.^2 + dy .^2));

vel_est = distances./diff(sample_times); % m/s
%disp(vel_est);

tot_dist = sum(distances);
%disp(tot_dist)

figure('Name','speed estimates')
plot(ty,v_peak, 'DisplayName','est') % km/h
hold on
plot(sample_times(1:end-1)+4, vel_est*3.6, 'DisplayName','GPS')
legend
grid on
ylabel('speed [km/h]')
xlabel('time [s]')

%%

fig = gcf;
exportgraphics(fig, 'GPSMap_01.pdf', 'ContentType', 'vector');


%%


Gx_trunc = Gx(wsize/2:end-wsize/2);
Gy_trunc = Gy(wsize/2:end-wsize/2);

thet = tht(:);
    
heading_rate = Gx_trunc.*cos(thet)+Gy_trunc.*sin(thet);
roll_rate    = -Gx_trunc.*sin(thet)+Gy_trunc.*cos(thet);

figure("Name","heading_rate")
plot(ty,heading_rate);
title('heading rate');

fc = 1; 
fs = 100;
n = 100; % filter order
b = fir1(n, (fc/(fs/2)), 'low');

heading_filtered = filtfilt(b,1,heading_rate);

vx = v_peak.*cos(heading_filtered);
vy = v_peak.*sin(heading_filtered);

x = cumtrapz(ty, vx);
y = cumtrapz(ty, vy);

figure;
plot(x,y); axis equal
title('x, y')

%%

X = [cos(thet), sin(thet), ones(size(thet))];  % model basis
coeff_x = X \ Gx_trunc;    % [Ax; Bx; offset_x]
coeff_y = X \ Gy_trunc;    % [Ay; By; offset_y]

Ax = coeff_x(1); Bx = coeff_x(2); ox = coeff_x(3);
Ay = coeff_y(1); By = coeff_y(2); oy = coeff_y(3);

amp_x = sqrt(Ax^2 + Bx^2);
amp_y = sqrt(Ay^2 + By^2);

fprintf('Leakage amplitudes: Gx: %.4f, Gy: %.4f\n', amp_x, amp_y);