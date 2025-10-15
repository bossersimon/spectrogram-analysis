
%% latest estimation method

clear all
accelScale = 1;%/9.82;

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

x = M_filt(:,1);
y = M_filt(:,2);

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
[sx,fx,tx, Px] = spectrogram(x,win,ovlap,Ndft,fs);
[sy,fy,ty, Py] = spectrogram(y,win,ovlap,Ndft,fs);


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


% %% Parameter estimation using STFT 

% f_vals = fy(f0_idx); % values of largest peaks in the sampled spectrum
% 
% cols = 1:size(sy,2);   % cols for sub2ind
% Fky = sy(sub2ind(size(sy), f0_idx, cols))/sum(win);  % 
% Fkx = sx(sub2ind(size(sx), f0_idx, cols))/sum(win);
% 
% %figure; hold on;
% %plot(ty, real(Fky) ) % real-valued DFT
% %plot(ty, real(Fkx) )
% %plot(real(Fkx),real(Fky))
% %title("real(Fky)")

%%
dt = 1/100;
N = size(M,1);   

% NEW threshold logic
gyro_vals = -Gz(wsize/2:end-wsize/2)/360; % DPS/360 [1/s]
freq_res = fs/Ndft;

f_th = 2; % Roughly 6 km/h
th_idx = round(f_th / freq_res) + 1;

f0_idx = nan(1,size(sy,2));

%%

curr_max_idx_y = 1;
curr_max_idx_x = 1;

bp_width = 0.3; % in Hz
max_jump = 5;
for t2 = 1:size(sy,2)
    if gyro_vals(t2)<f_th 
        % use gyro derived freq
        idx = round(gyro_vals(t2) / freq_res) + 1;
        idx = max(1, min(idx, Ndft));

        if abs(idx - curr_max_idx_y) <= max_jump
            f0_idx(t2) = idx;
            curr_max_idx_y = idx;
            curr_max_idx_x = idx;
        else
            f0_idx(t2) = curr_max_idx_y;  % ignore jump
        end
    else
        % obtain bp range from prev. window
        curr_max_idx = round((curr_max_idx_y + curr_max_idx_x)/2)+1;

        lo_f = max(0, fy(curr_max_idx) - bp_width);
        hi_f = fy(curr_max_idx) + bp_width;
    
        f_pass = fy >= lo_f & fy <= hi_f;
        sy_bp = zeros(size(sy(:,t2)));
        sx_bp = zeros(size(sy(:,t2)));

        sy_bp(f_pass) = sy(f_pass, t2);
        sx_bp(f_pass) = sx(f_pass, t2);
                    
        % average peak of x and y
        [~,max_idx_x] = max(abs(sx_bp));
        [~,max_idx_y] = max(abs(sy_bp));

        f0_idx(t2) = round((max_idx_x + max_idx_y)/2); % store
        curr_max_idx_y = max_idx_y; % update
        curr_max_idx_x = max_idx_x;
    end
end

f_vals = fy(f0_idx); % values of largest peaks in the sampled spectrum

%%
figure;
plot(ty, f0_idx, 'DisplayName','f0_idx')

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
colormap(bone);
colorbar;

colors = {'k', [0.3 0.3 0.3], [0.6 0.6 0.6], [0.9 0.9 0.9]};
styles = {'-', '--', '-.', ':'};

hold on;
plot(ty, v_lo, 'Color', colors{1},'LineStyle', styles{1}, 'LineWidth', 1.5);  % Lower band edge
plot(ty, v_hi, 'Color', colors{1}, 'LineWidth', 1.5);  % Upper band edge
plot(ty, v_peak, 'r-', 'LineWidth', 2);   % Tracked speed

gz = -Gz(wsize/2:end-wsize/2)*wheel_circ/100; % (DPS/360)*circ*3.6 [km/h]
p1 = plot(ty,gz,'Color',[1.0, 0.4, 0.0]);

% plot(ty, v_lo2, 'w--', 'LineWidth', 1.5);  % Lower band edge
% plot(ty, v_hi2, 'w--', 'LineWidth', 1.5);  % Upper band edge
% plot(ty, v_peak2, 'r-', 'LineWidth', 2);   % Tracked speed
legend('Lower Band Edge', 'Upper Band Edge', 'Tracked Speed','FontSize',12);

%%

% Peak idx stored in f0_idx

f_vals = fy(f0_idx);

frame_idx = 4020;
%frame_idx = 8506;

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
plot(t0, win.*x(win_idx),'DisplayName',"x")
hold on
plot(t0, win.*y(win_idx),'DisplayName',"x")
%plot(ty(frame_idx-wsize/2:frame_idx+wsize/2-1), win.*Ay(frame_idx-wsize/2:frame_idx+wsize/2-1), "DisplayName", 'y')
legend
grid on

%%
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

fq = (f0_idx(frame_idx)-1 + xq)*fs/Ndft; % convert to frequency

% Plot
plot(fq, yqy, 'r--', 'LineWidth', 1.5)  % quadratic interpolation curve
plot(fq, yqx, 'r--', 'LineWidth', 1.5)

yq_avg = (yqy+yqx)/2;
plot(fq, yq_avg, 'r--', 'LineWidth', 1.5)

%%

Nframes = length(f_interpy);
Fy_interp = zeros(1,Nframes);
Fx_interp = zeros(1,Nframes);

L = length(win);
n = (0:(L-1)).';

for m = 1:Nframes
    start_idx = m;
    idx = start_idx:(start_idx + wsize-1);

    y_frame = y(idx).*win(:);
    x_frame = x(idx).*win(:);

    f0 = f_interpy(m);

    Fy_interp(m) = y_frame' *exp(-1j*2*pi*f0.*n/fs);
    Fx_interp(m) = x_frame' *exp(-1j*2*pi*f0.*n/fs);
end

S = Fx_interp + 1j*Fy_interp;
tol = 1e-6;    % remove tiny noise
S(abs(S) < tol) = 0;
phi = angle(S);

phi_offs = transpose(2*pi*f_interpy*wsize/(2*fs));
theta = unwrap(phi)+phi_offs;
tht = wrapToPi(theta);

%%
plot(f_interpy)
hold on
plot(f_vals)

%% Diff. angle
%delta_phi =  angle(Fy_interp(2:end) .* conj(Fy_interp(1:end-1)) );
figure;
t_phi = t(wsize/2:end-wsize/2);
plot(t_phi,tht)

%%
%plot(t_phi, unwrap(tht))
plot(t_phi(1:end-1), diff(tht))
%%
hops = 1; 

delta_phi = unwrap(tht(hops+1:end)) - unwrap(tht(1:end-hops));

inst_freq = fs/(2*pi*hops) *delta_phi;
ty_hop = ty(1:end-hops);

figure; hold on;

plot(ty_hop, inst_freq, 'DisplayName', 'differencing', 'LineWidth',2)
plot(ty(1:end-1), f_interpy(1:end-1), 'DisplayName', 'interpolation','LineWidth',2)
plot(ty(1:end-1), f_vals(1:end-1) , 'DisplayName', 'no interpolation', 'LineWidth',2)
%plot(t,x/10+5)
%plot(t, x/10);
%plot(t,y/10);
p1 = plot(ty,gyro_vals,'DisplayName','Gyro');

title('Frequency estimates')
xlabel('Time [s]')
ylabel('Frequency [Hz]')
legend;
grid on

%% save figure
fig = gcf;
exportgraphics(fig, 'differencing03.pdf', 'ContentType', 'vector');

%%
%su = cumsum(delta_phi);
plot(t_phi ,unwrap(phi), 'DisplayName','dphi sum')
hold on
Fkx_angle = angle(Fkx);
plot(ty, unwrap(Fkx_angle) ,'DisplayName', 'unwr angle(Fxk)') 
legend


%%

% figure;
% plot(real(Fk), imag(Fk), '-o');
% axis equal;
% xlabel('Real Part');
% ylabel('Imaginary Part');
% title('Complex Trajectory of DFT Coefficient Over Time');

n = length(Fkx);
k_step = 1;  % color update interval (e.g., every 5 frames)
intervals = 1:k_step:n;
if intervals(end) < n
    intervals = [intervals, n];  % include final point
end
num_segments = 500;
cmap = turbo(5);
figure; hold on;

for i = 1:num_segments
    k1 = intervals(i);
    k2 = intervals(i+1);
    plot(real(S(k1:k2)), imag(S(k1:k2)), '-','Color', cmap(mod(i,5)+1,:), 'LineWidth', 1.5);
end
plot(real(S(1)), imag(S(1)), 'ro') 
plot(real(S(k2)), imag(S(k2)), 'bo') 


grid on

axis equal;
xlabel('Real [rad/s]');
ylabel('Imag [rad/s]');
colormap(cmap);
