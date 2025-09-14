
%% lsq.m


accelScale = 1/9.82; % scale accelerometer readings


M = readmatrix("recordings/recording_20250701_02.csv");

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
tiledlayout(3,1,"TileIndexing","columnmajor")
ax = [];
for i = 1:3
    ax(end+1) =nexttile;
    if i<4
        plot(t, M(:,i)/accelScale)
        hold on
        plot(t,M_filt(:,i)/accelScale)
    else
        plot(t, M(:,i))
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

Ax = M_filt(:,1);
Ay = M_filt(:,2);
 
wheel_circ = 1.82;


%% Parameter estimation using lsq 

model = @(b,t) cos(2*pi*b(1)*t+b(2)) + b(3);
b0 = [0.1, 0.1, 0.1]; % initial guess

lb = [0, -pi, min(Ay)];
ub = [5.5, pi, max(Ay)];

opts = optimoptions('lsqcurvefit', ...
    'Algorithm','levenberg-marquardt', ...
    'MaxFunctionEvaluations', 5000, ...
     'Display', 'off'); % 'iter'

window_size = 10;
N = length(Ay);
step = 5;

% preallocation
n = floor((N - window_size) / step) + 1;
time = zeros(1,n);
xvals = zeros(1,n);
yvals = zeros(1,n);
xparams = cell(1,n);
yparams = cell(1,n);
t_windows = zeros(n,window_size);
x_fits = zeros(n,window_size);
y_fits = zeros(n,window_size);

k=1;
center_idx = window_size/2;

for i=window_size/2:step:N-window_size/2
    j = i-window_size/2+1:(i+window_size/2);
    t_win = t(j);
    Ay_win = Ay(j);
    Ax_win = Ax(j);

    %t_rel = t_win-t_win(center_idx);
    t_rel = t_win;
  
    %b_hatx = lsqcurvefit(model, b0, t_rel, Ax_win, lb, ub, opts);
    %b_haty = lsqcurvefit(model, b0, t_rel, Ay_win, lb, ub, opts);
    b_hatx = lsqcurvefit(model, b0, t_rel, Ax_win, [], [], opts);
    b_haty = lsqcurvefit(model, b0, t_rel, Ay_win, [], [], opts);

    y_fit = model(b_haty,t_rel);
    x_fit = model(b_hatx,t_rel);

    idx = center_idx;

    time(k) = t_win(idx);
    yvals(k) = y_fit(idx);
    xvals(k) = x_fit(idx);
    xparams{k} = b_hatx;
    yparams{k} = b_haty;

    t_windows(k,:) = t_win;
    x_fits(k,:) = x_fit;
    y_fits(k,:) = y_fit;

    k=k+1;
end


%%

num_pts = 5;
half_width = floor(num_pts/2);

start_idx = center_idx - half_width;
end_idx = center_idx + half_width;

n = size(t_windows, 1);

t_plot = zeros(n, num_pts);
x_plot = zeros(n, num_pts);
y_plot = zeros(n, num_pts);

figure;

for k=1:n
    t_plot(k,:) = t_windows(k,start_idx:end_idx);
    x_plot(k,:) = x_fits(k,start_idx:end_idx);
    y_plot(k,:) = y_fits(k,start_idx:end_idx);
    
    plot(t_plot(k,:), x_plot(k,:))
    hold on
    %plot(t_plot(k,:), x_plot(k,:),'ro')
end
plot(t,Ax)


%%

t_vec = reshape(t_plot.',1,[]);
x_vec = reshape(x_plot.',1,[]);
y_vec = reshape(y_plot.',1,[]);

plot(t_vec, x_vec)
hold on
plot(t_vec, y_vec)



%% Plot full signals

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

xoffs = xparam_matrix(:,3);
yoffs = yparam_matrix(:,3);

xcorr = xvals - xoffs.';
ycorr = yvals - yoffs.';

figure;
plot(time,xvals,'DisplayName','xhat')
hold on
plot(time,xcorr,'DisplayName','xhat_ nooffset')
%plot(time,yvals)
legend

figure;
plot(xvals,yvals)
hold on
plot(xcorr,ycorr)
%%

arg_x = xparam_matrix(:,2);
arg_y = yparam_matrix(:,2);

phase_raw = atan2(yvals,xvals); %
phase_corr = atan2(ycorr,xcorr);
phase_unwr = unwrap(phase_raw);


%% Plot 3 frames

for i =40:220
    plot(t_windows(i,:),x_fits(i,:))
    hold on
    plot(t_windows(i,:),y_fits(i,:))
   % plot(t_windows(i,window_size/2), x_fits(i,window_size/2),'o')
   % plot(t_windows(i,window_size/2), y_fits(i,window_size/2),'o')
end

    plot(time(40:220),ycorr(40:220))
    plot(time(40:220),xcorr(40:220))
   % plot(time,xcorr)


%%


figure;
plot(time,phase_corr,'DisplayName','phase (offset removed)')
hold on
plot(time,phase_raw,'DisplayName','raw phase')
%plot(time,arg_x,'DisplayName','omega*t+phi')
ylabel('Angle [rad]')
xlabel('Time [s]')
title('LSQ Angle Estimate')
legend

fprintf('LSQ distance estimate: %.2f\n', phase_unwr(end)*wheel_circ/(2*pi));