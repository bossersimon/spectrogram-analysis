
% fig = figure('Units','pixels','Position',[100 100 800 600]);

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
by = fir1(n, (fc/(fs/2)), 'low');

M_filt = M;
M_filt(:,1:3) = filtfilt(by,1,M(:,1:3));

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
wsize = 100;
ovlap = wsize-1;
Ndft = 1024;

win = hann(wsize); % window function can be changed to something else
%win = gausswin(wsize);
%win = rectwin(wsize);
tic
[sx,fx,tx, Px] = spectrogram(Ax,win,ovlap,Ndft,fs);
[sy,fy,ty, Py] = spectrogram(Ay,win,ovlap,Ndft,fs);
toc

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
cy = colorbar;
cy.Label.String = 'Power/frequency [dB/Hz]';

hold on;
gyro_vals = -Gz(wsize/2:end-wsize/2)*wheel_circ/100; % (DPS/360)*circ*3.6 [km/h]
p1 = plot(ty,gyro_vals,'Color',[1.0, 0.4, 0.0]);
legend(p1, 'Gyroscope signal overlay', 'Location', 'northwest','FontSize', 12);


%%

PdB = 10*log10(Py);

sz = size(PdB,2);
indices = sz-700:sz;
f_indices = 1:150;

%%
% --- make figure ---
fig = figure('Units','normalized','OuterPosition',[0 0 1 1]);
imagesc(ty(indices), fy(f_indices), PdB(f_indices,indices));  % fx in kHz if desired
axis xy;                   % ensure correct orientation
colormap(jet);
axis tight;                % remove extra margins
axis off;                  % turn off all axes, ticks, etc.
title('');                 % no title
colorbar off;              % remove colorbar (we'll handle in pgfplots)

hold on
%gyro_vals = -Gz(wsize/2:end-wsize/2)/360; % (DPS/360)*circ*3.6 [km/h]
%p1 = plot(ty,gyro_vals,'Color',[1.0, 0.4, 0.0]);
%legend(p1, 'Gyroscope signal overlay', 'Location', 'northwest','FontSize', 12);

% --- export vector graphic ---
fname = 'myspectrogram_with_overlay';
print('-depsc2', fname);                 % cropped EPS export
system(['convert ' fname '.eps ' fname '.pdf']);  % use ImageMagick

% --- print axis and color limits ---
ax = gca;
fprintf('xmin: %g\nxmax: %g\n', ax.XLim);
fprintf('ymin: %g\nymax: %g\n', ax.YLim);
fprintf('zmin: %g\nzmax: %g\n', ax.CLim);


%% tikz plot

gyro_vals = -Gz(wsize/2:end-wsize/2)/360; % (DPS/360)*circ*3.6 [km/h]
p1 = plot(ty(indices),gyro_vals(indices),'Color',[1.0, 0.4, 0.0]);


%% Parameter estimation using STFT 

dt = 1/100;
N = size(M,1);   

% NEW threshold logic
gyro_vals = -Gz(wsize/2:end-wsize/2)/360; % DPS/360 [1/s]
freq_res = fs/Ndft;

f_th = 0; % Roughly 6 km/h
th_idx = round(f_th / freq_res) + 1;

f0_idx = nan(1,size(sy,2));

%%

curr_max_idx_y = 1;
curr_max_idx_x = 1;

bp_width = 0.3; % in Hz
max_jump = 5;
for t = 1:size(sy,2)
    if gyro_vals(t)<f_th 
        % use gyro derived freq
        idx = round(gyro_vals(t) / freq_res);
        idx = max(1, min(idx, Ndft));

        if abs(idx - curr_max_idx_y) <= max_jump
            f0_idx(t) = idx;
            curr_max_idx_y = idx;
            curr_max_idx_x = idx;
        else
            f0_idx(t) = curr_max_idx_y;  % ignore jump
        end
    else
        % obtain bp range from prev. window
        curr_max_idx = round((curr_max_idx_y + curr_max_idx_x)/2)+1;

        lo_f = max(0, fy(curr_max_idx) - bp_width);
        hi_f = fy(curr_max_idx) + bp_width;
    
        f_pass = fy >= lo_f & fy <= hi_f;
        sy_bp = zeros(size(sy(:,t)));
        sx_bp = zeros(size(sy(:,t)));

        sy_bp(f_pass) = sy(f_pass, t);
        sx_bp(f_pass) = sx(f_pass, t);
                    
        % average peak of x and y
        [~,max_idx_x] = max(abs(sx_bp));
        [~,max_idx_y] = max(abs(sy_bp));

        f0_idx(t) = round((max_idx_x + max_idx_y)/2); % store
        curr_max_idx_y = max_idx_y; % update
        curr_max_idx_x = max_idx_x;
    end
end

%%
figure;
plot(ty, f0_idx, 'DisplayName','f0_idx')
% hold on 
% plot(ty, f0_idx_dyn, 'DisplayName','f0_idx_dyn')
% legend
%%

f_peak = fy(f0_idx);        % Hz
v_peak = f_peak * wheel_circ * 3.6;

v_lo = max(v_peak - bp_width*wheel_circ*3.6, 0); 
v_hi = v_peak + bp_width*wheel_circ*3.6;

fig = figure('Units','normalized','OuterPosition',[0 0 1 1]); 
set(fig, 'PaperOrientation', 'landscape');

imagesc(ty, v_car, 10*log10(Py));
axis xy;
xlabel('Time (s)');
ylabel('Speed (km/h)');
title('Spectrogram');
%colormap(bone);
colormap(jet);
colorbar;

colors = {'k', [0.3 0.3 0.3], [0.6 0.6 0.6], [0.9 0.9 0.9]};
styles = {'-', '--', '-.', ':'};

hold on;
plot(ty, v_lo, 'Color', colors{1},'LineStyle', styles{1}, 'LineWidth', 1.5);  % Lower band edge
plot(ty, v_peak, 'r-', 'LineWidth', 2);   % Tracked speed
gyro_vals = -Gz(wsize/2:end-wsize/2)*wheel_circ/100; % (DPS/360)*circ*3.6 [km/h]
p1 = plot(ty,gyro_vals,'Color',[1.0, 0.4, 0.0]);
plot(ty, v_hi, 'Color', colors{1}, 'LineWidth', 1.5);  % Upper band edge

% plot(ty, v_lo2, 'w--', 'LineWidth', 1.5);  % Lower band edge
% plot(ty, v_hi2, 'w--', 'LineWidth', 1.5);  % Upper band edge
% plot(ty, v_peak2, 'r-', 'LineWidth', 2);   % Tracked speed
legend('Band Edges', 'Tracked Speed', 'Gyro', 'FontSize',12);

%%

indices = sz-600:sz;
f_indices = 1:150;

fig = figure('Units','normalized','OuterPosition',[0 0 1 1]);
imagesc(ty(indices), fy(f_indices), PdB(f_indices,indices));  % fx in kHz if desired
axis xy;                   % ensure correct orientation
colormap(jet);
axis tight;                % remove extra margins
axis off;                  % turn off all axes, ticks, etc.
title('');                 % no title
colorbar off;              % remove colorbar (we'll handle in pgfplots)

hold on
%gyro_vals = -Gz(wsize/2:end-wsize/2)/360; % (DPS/360)*circ*3.6 [km/h]
%p1 = plot(ty,gyro_vals,'Color',[1.0, 0.4, 0.0]);
%legend(p1, 'Gyroscope signal overlay', 'Location', 'northwest','FontSize', 12);

% --- export vector graphic ---
fname = 'myspectrogram2';
print('-depsc2', fname);                 % cropped EPS export
system(['convert ' fname '.eps ' fname '.pdf']);  % use ImageMagick

% --- print axis and color limits ---
ax = gca;
fprintf('xmin: %g\nxmax: %g\n', ax.XLim);
fprintf('ymin: %g\nymax: %g\n', ax.YLim);
fprintf('zmin: %g\nzmax: %g\n', ax.CLim);



%%
f_peak = fy(f0_idx);        % Hz
%v_peak = f_peak * wheel_circ * 3.6;

v_lo = max(f_peak - bp_width, 0); 
v_hi = f_peak + bp_width;

colors = {'k', [0.3 0.3 0.3], [0.6 0.6 0.6], [0.9 0.9 0.9]};
styles = {'-', '--', '-.', ':'};

hold on;
plot(ty(indices), v_lo(indices), 'Color', colors{1},'LineStyle', styles{1}, 'LineWidth', 1.5);  % Lower band edge
plot(ty(indices), f_peak(indices), 'r-', 'LineWidth', 2);   % Tracked speed
gyro_vals = -Gz(wsize/2:end-wsize/2)/360; % (DPS/360)*circ*3.6 [km/h]
p1 = plot(ty(indices),gyro_vals(indices),'Color',[1.0, 0.4, 0.0]);
plot(ty(indices), v_hi(indices), 'Color', colors{1}, 'LineWidth', 1.5);  % Upper band edge

legend('Band Edges', 'Tracked Speed', 'Gyro', 'FontSize',12);

%%
% Peak idx stored in f0_idx

f_vals = fy(f0_idx);

%frame_idx = 250;
frame_idx = 4200;

Fky_frame = sy(:, frame_idx)/sum(win);
Fkx_frame = sx(:, frame_idx)/sum(win);

subplot(2,1,1);
plot(fy(1:100), abs(Fky_frame(1:100)), 'b', "DisplayName", "y")
hold on
plot(fx(1:100), abs(Fkx_frame(1:100)), 'r', "DisplayName","x")
legend
grid on

%plot(f0_idx(frame_idx), f_vals(frame_idx), 'r*') 
plot(f_vals(frame_idx), abs(Fky_frame(f0_idx(frame_idx))), 'b*');
plot(f_vals(frame_idx), abs(Fkx_frame(f0_idx(frame_idx))), 'r*');

win_idx = frame_idx : (frame_idx + wsize - 1);

t0 = (win_idx - 1) / fs;
%t0 = (ty(frame_idx) - (wsize/2)/fs) : (1/fs) : (ty(frame_idx) + (wsize/2 - 1)/fs);
subplot(2,1,2);
%plot(tx(frame_idx-wsize/2:frame_idx+wsize/2-1), win.*Ax(frame_idx-wsize/2:frame_idx+wsize/2-1),'DisplayName',"x")
plot(t0, win.*Ax(win_idx),'DisplayName',"x")
hold on
plot(t0, win.*Ay(win_idx),'DisplayName',"x")
%plot(ty(frame_idx-wsize/2:frame_idx+wsize/2-1), win.*Ay(frame_idx-wsize/2:frame_idx+wsize/2-1), "DisplayName", 'y')
legend
grid on


%% Quadratic spectral peak interpolation

cols = 1:size(sy, 2);

alphay = abs(sy(sub2ind(size(sy), max(1,f0_idx-1), cols)))/sum(win); 
betay  = abs(sy(sub2ind(size(sy), f0_idx, cols)))/sum(win);  
gammay = abs(sy(sub2ind(size(sy), f0_idx+1, cols)))/sum(win); 

alphax = abs(sx(sub2ind(size(sy), max(1,f0_idx-1), cols)))/sum(win); 
betax  = abs(sx(sub2ind(size(sy), f0_idx, cols)))/sum(win);  
gammax = abs(sx(sub2ind(size(sy), f0_idx+1, cols)))/sum(win);

py = 0.5*(alphay-gammay)./(alphay- 2*betay + gammay); % fractional offset
px = 0.5*(alphax-gammax)./(alphax- 2*betax + gammax);

k_interpy = f0_idx(:) + py(:); % peak location in fractional bins
k_interpx = f0_idx(:) + px(:);

f_interpy = max((k_interpy(:)-1)*fs/Ndft,0); % corresponding peak frequency
f_interpx = max((k_interpx(:)-1)*fs/Ndft,0);
f_interp_avg = (f_interpy+f_interpx)/2;

max_interpy = betay - (1/4) * (alphay - gammay) .* py; % peak magnitude estimate
max_interpx = betax - (1/4) * (alphax - gammax) .* px;
max_interp_avg = (max_interpx+max_interpy)/2;

plot(fy(1:100), abs(Fky_frame(1:100)))
hold on
plot(fy(1:100), abs(Fkx_frame(1:100)), 'r')

plot(fy(f0_idx(frame_idx)), abs(Fky_frame(f0_idx(frame_idx))), 'go'); % bin peak
plot(f_interpy(frame_idx), max_interpy(frame_idx) , 'bo') % <---interpolated peak
plot(f_interpx(frame_idx), max_interpx(frame_idx) , 'bo')
plot(f_interp_avg(frame_idx), max_interp_avg(frame_idx), 'pentagram', 'Color', 'm')

plot(fy(f0_idx(frame_idx)-1), alphay(frame_idx), 'r*')
%plot(fy(f0_idx(frame_idx)), beta(frame_idx), 'r*')
plot(fy(f0_idx(frame_idx)+1), gammay(frame_idx), 'r*')

plot(fx(f0_idx(frame_idx)-1), alphax(frame_idx), 'r*')
%plot(fy(f0_idx(frame_idx)), beta(frame_idx), 'r*')
plot(fx(f0_idx(frame_idx)+1), gammax(frame_idx), 'r*')


% parabola coefficients
ay = 0.5*(alphay(frame_idx) - 2*betay(frame_idx) + gammay(frame_idx));
by = 0.5*(gammay(frame_idx) - alphay(frame_idx));
cy = betay(frame_idx);

ax = 0.5*(alphax(frame_idx) - 2*betax(frame_idx) + gammax(frame_idx));
bx = 0.5*(gammax(frame_idx) - alphax(frame_idx));
cx = betax(frame_idx);

xq = linspace(-1.2, 1.2, 200);  % fine resolution around the peak
yqy = ay*xq.^2 + by*xq + cy;        % parabola
yqx = ax*xq.^2 + bx*xq + cx; 

frame_idx = 4200;


fq = (f0_idx(frame_idx)-1 + xq)*fs/Ndft; % convert to frequency

% Plot
plot(fq, yqy, 'r--', 'LineWidth', 1.5)  % quadratic interpolation curve
plot(fq, yqx, 'r--', 'LineWidth', 1.5)

yq_avg = (yqy+yqx)/2;
plot(fq, yq_avg, 'r--', 'LineWidth', 1.5)


%% Obtaining the phase at interpolated peak

Nframes = length(f_interpy);
Fy_interp = zeros(1,Nframes);
Fx_interp = zeros(1,Nframes);

L = length(win);
n = (0:(L-1)).';

for m = 1:Nframes
    start_idx = m;
    idx = start_idx:(start_idx + wsize-1);

    y_frame = Ay(idx).*win(:);
    x_frame = Ax(idx).*win(:);

    f0 = f_interpy(m);

    Fy_interp(m) = y_frame' *exp(-1j*2*pi*f0.*n/fs);
    Fx_interp(m) = x_frame' *exp(-1j*2*pi*f0.*n/fs);
end
%%

model_dft = @(t,b) cos(b{1}) + b{2};

% cols = 1:size(sy,2);   % cols for sub2ind
% %Sx = sx(sub2ind(size(sx), f0_idx, cols));  % values of S at wanted frequencies
% %Sy = sy(sub2ind(size(sy), f0_idx, cols));
% 
% % phase offset estimate (phi_0)
% S = Sx + 1j*Sy;
S = Fx_interp + 1j*Fy_interp;
% 
tol = 1e-6;    % remove tiny noise
S(abs(S) < tol) = 0;
% 
% %theta_dft = unwrap(angle(S));
theta_dft = angle(S);

phi_offs = transpose(2*pi*f_interpy*wsize/(2*fs));

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
%xlabel(tl,'Time [s]','FontSize', 20, 'Interpreter', 'latex');
% 
% if accelScale < 1
%     ylabel(tl,"Acceleration [m/s$^2$]", 'Interpreter','latex', 'FontSize', 18);
% else
%     ylabel(tl,"Acceleration [g]", 'Interpreter','latex', 'FontSize', 18);
% end

ax(end+1) =nexttile;
plot(t, Ax/accelScale,'DisplayName','raw x','LineWidth',1.5)
hold on
plot(tx,xhat/accelScale,'DisplayName','xhat','LineWidth',1.5)
ylabel(ax,"Acceleration [m/s$^2$]", 'Interpreter','latex', 'FontSize', 18);

title('x-axis','FontSize', 20, 'Interpreter', 'latex')
grid on
legend('FontSize', 16)

ax(end+i)=nexttile;
plot(t, Ay/accelScale,'DisplayName','raw y','LineWidth',1.5)
hold on
plot(ty,yhat/accelScale,'DisplayName','yhat','LineWidth',1.5)
ylabel(ax(end),"Acceleration [m/s$^2$]", 'Interpreter','latex', 'FontSize', 18);
xlabel(ax(end),'Time [s]','FontSize', 18, 'Interpreter', 'latex');

title('y-axis','FontSize', 20, 'Interpreter', 'latex')
grid on
legend('FontSize', 16)

linkaxes(ax,"x")

%%
fig = gcf;
exportgraphics(fig, 'spectrogram2.pdf', 'ContentType', 'vector');

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
%plot(ty, tht,'mo','DisplayName','phase')
ylabel("Phase [rad]", 'FontSize', 18, 'Interpreter','latex');
title('Phase', 'FontSize', 20, 'Interpreter', 'latex');
grid on

linkaxes(ax,"x")
%% Total distance:

plot(ty,tht)
d = unwrap(tht);
d = d(end)*wheel_circ/(2*pi);

disp([d]);

%% export

fig = gcf;
exportgraphics(fig, 'phasedist_02.pdf', 'ContentType', 'vector');


%% 

fc = 0.5; 
fs = 100;
n = 100; % filter order
by = fir1(n, (fc/(fs/2)), 'high');

x_filt = filtfilt(by,1,Ax);
y_filt = filtfilt(by,1,Ay);



figure;
plot(xhat, yhat)
hold on
%plot(Ax,Ay)

%%
theta_raw  = atan2(Ay,Ax);
theta_raw2  = atan2(y_filt,x_filt);

fig = figure('Units','normalized','OuterPosition',[0 0 1 1]); 
set(fig, 'PaperOrientation', 'landscape');

%plot(ty,theta_raw, 'DisplayName','raw no offs')
hold on
%plot(ty,tht, 'DisplayName','est')
plot(t,theta_raw, 'DisplayName','raw','Color','b', 'LineWidth',1.5)
title("raw phase")
xlabel('Time [s]')
ylabel('Angle [rad]')
grid on
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


gpslog = readmatrix("recordings/GPS/sensorLog_20250701_02.txt");

time= gpslog(:,1);
lat = gpslog(:,3);
lon = gpslog(:,4);
alt = gpslog(:,5);


fig = figure('Units','normalized','OuterPosition',[0 0 1 1]); 
set(fig, 'PaperOrientation', 'landscape');

geoplot(lat,lon,"-.",'LineWidth',1,'Color','r')
geobasemap streets
title("Car route")

%%


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
%disp(sample_times);

distances = sqrt(dx.^2 + dy.^2);

vel_est = distances./diff(sample_times); % m/s
%disp(vel_est);

tot_dist = sum(distances);
disp("tot_dist: " + tot_dist)

%%
figure('Name','speed estimates')
plot(ty,v_peak, 'DisplayName','Tracked Speed','LineWidth',2) % km/h
hold on
plot(sample_times(1:end-1)-14, vel_est*3.6, 'DisplayName','GPS approximate','LineWidth',2)
legend
grid on
ylabel('speed [km/h]')
xlabel('time [s]')

%%

fig = gcf;
exportgraphics(fig, 'GPSspeed.pdf', 'ContentType', 'vector');


%% This does not seem to work

Gx_trunc = Gx(wsize/2:end-wsize/2);
Gy_trunc = Gy(wsize/2:end-wsize/2);

thet = tht(:);
    
heading_rate = Gx_trunc.*cos(thet) + Gy_trunc.*sin(thet);
%roll_rate    = -Gx_trunc.*sin(thet)+Gy_trunc.*cos(thet);

figure("Name","heading_rate")
plot(ty,heading_rate);
hold on
title('heading rate');

vx = v_peak.*cos(heading_rate);
vy = v_peak.*sin(heading_rate);

x = cumtrapz(ty, vx);
y = cumtrapz(ty, vy);

%figure;
%plot(x,y); axis equal
%title('x, y')