
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

%win = hann(wsize); % window function can be changed to something else
win = gausswin(wsize);
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


%% Parameter estimation using STFT 

dt = 1/100;
N = size(M,1);   

% NEW threshold logic
gyro_vals = -Gz(wsize/2:end-wsize/2)/360; % DPS/360 [1/s]
freq_res = fs/Ndft;

f_th = 1.7; % Roughly 6 km/h
th_idx = round(f_th / freq_res) + 1;

f0_idx = nan(1,size(sy,2));  % bin indices of largest peaks in the sampled spectrum

curr_max_idx = 1;

%Py_bp = zeros(size(sy));
bp_width = 0.3; % in Hz
max_jump = 3;
for t = 1:size(sy,2)
    if gyro_vals(t)<f_th 
        % center around gyro freq
        
        idx = round(gyro_vals(t) / freq_res) + 1;
        idx = max(1, min(idx, Ndft));
        
        if abs(idx - curr_max_idx) <= max_jump
            f0_idx(t) = idx;
            curr_max_idx = idx;
       % else
      %      f0_idx(t) = curr_max_idx;  % ignore jump
        end
    end
    %else
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
    %end
end

f_vals = fy(f0_idx); % values of largest peaks in the sampled spectrum

cols = 1:size(sy,2);   % cols for sub2ind
Fky = sy(sub2ind(size(sy), f0_idx, cols))/sum(win);  % 
Fkx = sx(sub2ind(size(sx), f0_idx, cols))/sum(win);

figure;
plot(ty, real(Fky) ) % real-valued DFT
%plot(real(Fkx),real(Fky))
title("real(Fky)")
hold on

%%

% Display one DFT frame.
% Peak idx stored in f0_idx
% Remaining issue is that the bin resolution is finite, and the peak
% frequency jumps between bins. 
% Should attempt to interpolate the spectrum.

frame_idx = 4020;

Fky_frame = sy(:, frame_idx) / sum(win);
Fkx_frame = sx(:, frame_idx) / sum(win);

plot(fy(1:100), abs(Fky_frame(1:100)))
hold on
%plot(f0_idx(frame_idx), f_vals(frame_idx), 'r*') 
plot(f_vals(frame_idx), abs(Fky_frame(f0_idx(frame_idx))), 'r*');



%% Quadratic spectral peak interpolation

cols = 1:size(sy, 2);
alpha = abs(sy(sub2ind(size(sy), max(1,f0_idx-1), cols)))/sum(win); 
beta  = abs(sy(sub2ind(size(sy), f0_idx, cols)))/sum(win);  
gamma = abs(sy(sub2ind(size(sy), f0_idx+1, cols)))/sum(win); 


p = 0.5*(alpha-gamma)./(alpha- 2*beta + gamma); % fractional offset 

k_interp = f0_idx(:) + p(:); % peak location in fractional bins

f_interp = (k_interp(:)-1)*fs/Ndft; % corresponding peak frequency

max_interp = beta - (1/4) * (alpha - gamma) .* p; % peak magnitude estimate

plot(fy(1:100), abs(Fky_frame(1:100)))
hold on
plot(fy(f0_idx(frame_idx)), abs(Fky_frame(f0_idx(frame_idx))), 'go');
plot(f_interp(frame_idx), max_interp(frame_idx) , 'bo') % <---interpolated peak

plot(fy(f0_idx(frame_idx)-1), alpha(frame_idx), 'r*')
%plot(fy(f0_idx(frame_idx)), beta(frame_idx), 'r*')
plot(fy(f0_idx(frame_idx)+1), gamma(frame_idx), 'r*')


% parabola coefficients
a = 0.5*(alpha(frame_idx) - 2*beta(frame_idx) + gamma(frame_idx));
b = 0.5*(gamma(frame_idx) - alpha(frame_idx));
c = beta(frame_idx);

xq = linspace(-1.2, 1.2, 200);  % fine resolution around the peak
yq = a*xq.^2 + b*xq + c;        % parabola

fq = (f0_idx(frame_idx)-1 + xq)*fs/Ndft; % convert to frequency

% Plot
plot(fq, yq, 'r--', 'LineWidth', 1.5)  % quadratic interpolation curve


%% Obtaining the phase at interpolated peak

Nframes = length(f_interp);
Fy_interp = zeros(1,Nframes);
L = length(win);
n = (0:(L-1)).'; 

for m = 1:Nframes
    start_idx = m;
    idx = start_idx:(start_idx + wsize-1);

    y_frame = Ay(idx).*win(:);
    f0 = f_interp(m);

    Fy_interp(m) = y_frame' *exp(-1j*2*pi*f0.*n/fs);
end

%% Diff. angle

%delta_phi = angle(Fy_interp(2:end) ./ Fy_interp(1:end-1));  % size: 1 x (N-1)
delta_phi =  angle(Fy_interp(2:end) .* conj(Fy_interp(1:end-1)) );

%other_angle = atan2(imag(Fky),real(Fkx));
inst_freq = fs/(2*pi) *delta_phi;
plot(ty(1:end-1), inst_freq, 'DisplayName', 'f_inst')
hold on
plot(ty(1:end-1), f_interp(1:end-1), 'DisplayName', 'f interpolated')
plot(ty(1:end-1), f_vals(1:end-1) , 'DisplayName', 'no interpolation')
plot(ty(1:end-1), angle(Fy_interp(1:end-1)), 'DisplayName', 'phase1')
plot(ty(2:end), angle(Fy_interp(1:end-1)), 'DisplayName', 'phase2')

title('Frequency estimates')
legend;

%%
figure;
plot(ty, real(Fy_interp))


%%
su = cumsum(delta_phi);
plot(ty(1:end-1), cumsum(delta_phi), 'DisplayName','dphi sum')
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
k_step = 5;  % color update interval (e.g., every 5 frames)
intervals = 1:k_step:n;
if intervals(end) < n
    intervals = [intervals, n];  % include final point
end
num_segments = 200;
cmap = parula(num_segments);
figure; hold on;
for i = 1:num_segments
    k1 = intervals(i);
    k2 = intervals(i+1);
    plot(real(Fky(k1:k2)), imag(Fky(k1:k2)), '-','Color', cmap(i,:), 'LineWidth', 0.5);
    %plot(real(Fy_interp(k1:k2)), imag(Fy_interp(k1:k2)), '-','Color', cmap(i,:), 'LineWidth', 0.5);
end

axis equal;
xlabel('Real [rad/s]');
ylabel('Imag [rad/s]');
colormap(cmap);
