
%% lsq.m

clear all

accelScale = 1/9.82; % scale accelerometer readings


M = readmatrix("recordings/recording_20250701_06.csv");

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

%Ax = M_filt(:,1);
%Ay = M_filt(:,2);
Ax = M(:,1);
Ay = M(:,2);

 
wheel_circ = 1.82;


%% Parameter estimation using lsq 

model = @(b,t) cos(2*pi*b(1)*t+b(2)) + b(3);
b0 = [0.1, 0.1, 0.1]; % initial guess
b0x = b0;
b0y = b0;

lb = [0, 0, min(Ay)];
ub = [5.5, 2*pi, max(Ay)];

opts = optimoptions('lsqcurvefit', ...
    'Algorithm','levenberg-marquardt', ...
    'MaxFunctionEvaluations', 5000, ...
     'Display', 'none'); % 'iter'

window_size = 80;
N = length(Ay);
step = 40;

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

tic
for i=window_size/2:step:N-window_size/2
    j = i-window_size/2+1:(i+window_size/2);
    t_win = t(j);
    Ay_win = Ay(j);
    Ax_win = Ax(j);

    t_rel = t_win;

    b_hatx = lsqcurvefit(model, b0x, t_rel, Ax_win, [], [], opts);
    b_haty = lsqcurvefit(model, b0y, t_rel, Ay_win, [], [], opts);

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

    b0x = b_hatx;
    b0y = b_haty;

    k=k+1;
end
toc

%%

num_pts = 39;
half_width = floor(num_pts/2);

start_idx = center_idx - half_width;
end_idx = center_idx + half_width;

n = size(t_windows, 1);

t_plot = zeros(n, num_pts); % contains all individual fits
x_plot = zeros(n, num_pts);
y_plot = zeros(n, num_pts);


fig = figure('Units','normalized','OuterPosition',[0 0 1 1]); 
set(fig, 'PaperOrientation', 'landscape');

plot(t,Ax, 'Color', [0.6, 0.6, 0.6], 'LineWidth',1)
hold on
for k=1:n
    t_plot(k,:) = t_windows(k,start_idx:end_idx);
    x_plot(k,:) = x_fits(k,start_idx:end_idx);
    y_plot(k,:) = y_fits(k,start_idx:end_idx);
    
    plot(t_plot(k,:), x_plot(k,:), 'LineWidth',2, 'LineStyle','-');
    %plot(t_plot(k,:), x_plot(k,:),'ro')
end
grid on

%%

fig = gcf;
%exportgraphics(fig, 'lsq05.pdf', 'ContentType', 'vector');
%%

t_vec = reshape(t_plot.',1,[]);
x_vec = reshape(x_plot.',1,[]);
y_vec = reshape(y_plot.',1,[]);

figure;
plot(t_vec, x_vec)
hold on
plot(t_vec, y_vec)

%% Plot full signals

figure;
plot(t,Ay/accelScale, 'Color','k')
hold on
plot(t_vec, y_vec/accelScale,'Color', 'r')

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
plot(t_vec,x_vec/accelScale,'Color', 'r')

if accelScale < 1
    ylabel('Acceleration [m/s²]')
else
    ylabel('Acceleration [g]')
end

xlabel('Time [s]')
title('X-axis')
grid on



%% DC removal

xparam_matrix = cell2mat(xparams');
yparam_matrix = cell2mat(yparams');

xoffs = xparam_matrix(:,3);
yoffs = yparam_matrix(:,3);

x_corrected = x_plot - xoffs;
y_corrected = y_plot - yoffs;

x_corrected = reshape(x_corrected.',1,[]);
y_corrected = reshape(y_corrected.',1,[]);

figure;
plot(t_vec,x_vec,'DisplayName','xhat', 'LineWidth',2)
hold on
plot(t_vec,x_corrected,'DisplayName','xhat_ nooffset','LineWidth',2)
legend
grid on

%%

figure;
plot(t_vec, x_corrected)
hold on
plot(t_vec, y_corrected)

%%

% figure;
% plot(x_vec,y_vec)
% hold on
% plot(x_corrected,y_corrected)
% legend('estimate', 'estimate_nooffset')

%% phase from curve fit  

%arg_x = xparam_matrix(:,2);
%arg_y = yparam_matrix(:,2);

fx = abs(xparam_matrix(:,1));
fy = abs(yparam_matrix(:,1));

plot(fx, 'Color', [0 0 1 0.5])
hold on
plot(fy, 'Color', [0 1 0 0.5] )
plot((fx+fy)/2, 'Color', 'r')
grid on

xlabel('Time [s]')
ylabel('Rotational Frequency [Hz]')

%% phase from atan2

phase_corrected = atan2(y_corrected,x_corrected);


figure;
plot(t_vec,phase_corrected, 'LineWidth', 2)

%% plots both curves and phase

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
plot(t_vec,y_corrected,'DisplayName','yhat', 'LineWidth',1.5)

hold on
plot(t_vec,x_corrected,'DisplayName','xhat','LineWidth',1.5)

if accelScale < 1
    ylabel("Acceleration [m/s$^2$]", 'Interpreter','latex', 'FontSize', 18);
else
    ylabel("Acceleration [g]", 'Interpreter','latex', 'FontSize', 18);
end

title('Time signals', 'FontSize', 20, 'Interpreter', 'latex');
grid on
legend

ax(end+i)=nexttile;
plot(t_vec, phase_corrected,'DisplayName','phase','LineWidth',1.5,'Color','b')
hold on
grid on

linkaxes(ax,'x')

%% Distance calculation

unwr = unwrap(phase_corrected);
d = unwr(end)*wheel_circ/(2*pi);

disp("tot_dist: " +  d);



%%

fc = 1; 
fs = 100;
n = 100; % filter order
by = fir1(n, (fc/(fs/2)), 'high');

t = (0:N-1)*dt;
t= transpose(t);

Ax_filtered = filtfilt(by,1,M_filt(:,1));
Ay_filtered = filtfilt(by,1,M_filt(:,2));

figure;
plot(t,atan2(Ay_filtered, Ax_filtered))
hold on
plot(t_vec, phase)
legend('raw', 'fit')
