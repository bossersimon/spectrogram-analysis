
%% latest estimation method


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
l = (0:N-1)*dt;
l= transpose(l);
tiledlayout(3,2,"TileIndexing","columnmajor")
ax = [];
for i = 1:6
    ax(end+1) =nexttile;
    if i<4
        plot(l, M(:,i)/accelScale)
        hold on
        plot(l,M_filt(:,i)/accelScale)
    else
        plot(l, M(:,i))
    end

    if i == 6
        plot(l,Gz_corrected)
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
Ndft = 8192;

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

% NEW threshold logic
gyro_vals = -Gz(wsize/2:end-wsize/2)/360; % DPS/360 [1/s]
freq_res = fs/Ndft;

f_th = 1.7; % Roughly 6 km/h
th_idx = round(f_th / freq_res) + 1;

f0_idx = nan(1,size(sy,2));
for l = 1:size(sy,2)
    if gyro_vals(l)<f_th
        % use gyro derived freq
        idx = round(gyro_vals(l) / freq_res) + 1;
        idx = max(1, min(idx, Ndft));
        f0_idx(l) = idx;
    else
        % use STFT peak
        [~,rel_idx] = max(abs(sy(th_idx:end,l)));
        idx = rel_idx + th_idx - 1;
        f0_idx(l) = idx;
    end
end

f_vals = fy(f0_idx);

cols = 1:size(sy,2);   % cols for sub2ind
%Fk = sy(sub2ind(size(sy), f0_idx, cols))/sum(win);
Fk = sx(sub2ind(size(sx), f0_idx, cols))/sum(win);


figure;
plot(t(wsize/2:end-wsize/2),real(Fk))

%%

delta_phi = unwrap(angle(Fk(2:end) ./ Fk(1:end-1)));  % size: 1 x (N-1)

% figure;
% plot(real(Fk), imag(Fk), '-o');
% axis equal;
% xlabel('Real Part');
% ylabel('Imaginary Part');
% title('Complex Trajectory of DFT Coefficient Over Time');

n = length(Fk);
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
    plot(real(Fk(k1:k2)), imag(Fk(k1:k2)), 'o','Color', cmap(i,:), 'LineWidth', 2);
end

axis equal;
xlabel('Real');
ylabel('Imag');
colormap(cmap);