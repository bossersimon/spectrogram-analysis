%% lsq2.m

clear all 


t = 0:1/100:40;
f = ones(1,length(t));
f = 2*f;
phi_offs = 0;
Ax = cos(2*pi*f.*t) + 0.1*rand(size(t)) + 0.02*t;
Ay = sin(2*pi*f.*t) + 0.1*rand(size(t)) - 0.02*t;

plot(t,Ax);
hold on
plot(t,Ay);

%% 
% Ax = chirp(t,0,40,5)+0.05*rand(size(t));
% Ay = chirp(t,0,40,5)+0.05*rand(size(t));
% plot(t,Ax)

%%

model = @(b,t) cos(2*pi*b(1)*t+b(2)) + b(3);
b0 = [0.1, 0.1, 0.1]; % initial guess
b0x = b0;
b0y = b0;

lb = [0, -pi+1, min(Ay)];
ub = [5.5, pi, max(Ay)];

opts = optimoptions('lsqcurvefit', ...
    'Algorithm','levenberg-marquardt', ...
    'MaxFunctionEvaluations', 5000, ...
     'Display', 'off'); % 'iter'

window_size = 16;
N = length(Ay);
step = 8;

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

    t_rel = t_win-t_win(center_idx);
    %t_rel = t_win;

    %b_hatx = lsqcurvefit(model, b0, t_rel, Ax_win, lb, ub, opts);
    %b_haty = lsqcurvefit(model, b0, t_rel, Ay_win, lb, ub, opts);
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

%%
num_pts = 9;
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
    
    plot(t_plot(k,:), x_plot(k,:), 'LineWidth',2)
    hold on
    %plot(t_plot(k,:), x_plot(k,:),'ro')
end
plot(t,Ax,'LineStyle','-.')
%%

t_vec = reshape(t_plot.',1,[]);
x_vec = reshape(x_plot.',1,[]);
y_vec = reshape(y_plot.',1,[]);

%plot(t_vec, x_vec)
%hold on
%plot(t_vec, y_vec)


%% DC removal

xparam_matrix = cell2mat(xparams');
yparam_matrix = cell2mat(yparams');

Bx = xparam_matrix(:,3);
By = yparam_matrix(:,3);

x_corr = x_plot-Bx;
y_corr = y_plot-By;

x_vec = reshape(x_corr.',1,[]);
y_vec = reshape(y_corr.',1,[]);

plot(t_vec, x_vec)
hold on
%plot(t_vec, y_vec)


%% First method: "interpolate phase from frequency"

freqx = xparam_matrix(:,1);
freqy = yparam_matrix(:,1);
argx = xparam_matrix(:,2);
argy = yparam_matrix(:,2);

t_diff = t_plot(1,:) - t_plot(1,ceil(end/2));
phisx = zeros(size(t_plot));
phisy= zeros(size(t_plot));

for k=1:size(t_plot,1)
    phisy(k,:) = argy(k)+freqy(k)*2*pi*t_diff;
    phisx(k,:) = argx(k)+freqx(k)*2*pi*t_diff;
end


%%
yphi_vec = reshape(phisy.',1,[]);

plot(t_vec,yphi_vec)
figure;
plot(t_vec,unwrap(yphi_vec))

%% Second method: phase from atan2
phase = atan2(y_vec,x_vec);
plot(t_vec, unwrap(phase))


%%


figure;

plot(t,Ay, 'Color','k')
hold on
plot(t_vec,y_vec,'Color', 'r')

legend('Ay', 'Reconstructed')


% 
% figure;
% 
% plot(t,Ax, 'Color','k')
% hold on
% plot(time,xvals,'Color', 'r')
%%

figure;
plot(t_vec,x_vec)
hold on
plot(t_vec, y_vec)
legend('xhat', 'yhat')


%%

fig = figure('Units','normalized','OuterPosition',[0 0 1 1]); 
set(fig, 'PaperOrientation', 'landscape');

tl = tiledlayout(2,1,"TileIndexing","columnmajor");
xlabel(tl,'Time [s]','Interpreter','latex', 'FontSize', 20);

ax = [];
ax(end+1) =nexttile;
plot(t_vec,y_vec,'DisplayName','yhat', 'LineWidth',1.5)

hold on
plot(t_vec,x_vec,'DisplayName','xhat','LineWidth',1.5)

title('Time signals', 'FontSize', 20, 'Interpreter', 'latex');
grid on
legend

ax(end+i)=nexttile;
tht = wrapToPi(phase);

plot(t_vec, tht,'DisplayName','phase','LineWidth',1.5,'Color','b')
hold on
%plot(ty, tht,'mo','DisplayName','phase')
ylabel("Phase [rad]", 'FontSize', 18, 'Interpreter','latex');
title('Phase', 'FontSize', 20, 'Interpreter', 'latex');
grid on

linkaxes(ax,"x")


%%


argx_unwr = unwrap(argx);
argy_unwr = unwrap(argy);

figure;
plot(time, argx)
%plot(time,argx_unwr)
hold on
plot(time, argy)
%plot(time,argy_unwr)



%% alt. method

window_size = 10;
N = length(Ay);
step = 1;

n = floor((N - window_size) / step) + 1; % number of fits
time = zeros(1,n);

bx = cell(1,n);
by = zeros(3,n);
w_best = zeros(1,n);
%w = 1e-4; % start frequency
k=1;
for i=window_size/2:step:N-window_size/2 %  i is center idx
    j = i-window_size/2+1:(i+window_size/2); % same indices but shifted half a window
    t_win = t(j);
    Ay_win = Ay(j);
    %Ax_win = Ax(j);
    %t_rel = t_win-t_win(1);
    t_rel = t_win;

    freqs = linspace(0.01, 5, 200); % min,max,num
    err = zeros(1, length(freqs));
    best_beta = zeros(3,1);
    best_w_idx = 1;
    
    min_err = Inf;
    % try different frequencies
    for idx = 1:length(freqs)
        w = 2*pi*freqs(idx);
        % Csin, Dcos +B
        X = [sin(w*t_rel(:)), cos(w*t_rel(:)), ones(length(t_rel),1)];
        [Q, R] = qr(X, 0);    

        beta = R \ (Q.' * Ay_win(:));
        %beta = X\Ay_win(:);
        res = Ay_win(:)- X*beta;
        err(idx) = sum(res.^2);

        if err(idx) < min_err
            min_err = err(idx);
            best_beta = beta;
            best_w_idx = idx;
        end
        
        %bx{k} = X\Ax_win;

    end
    w_best(k) = 2*pi*freqs(best_w_idx);
    by(:,k) = best_beta;

    k=k+1;
end

%%

Cy = by(1,:).';
Dy = by(2,:).';
By = by(3,:).';

%scale_y = 1/sqrt(Cy0.^2+Dy0.^2);
%Cy = Cy0.*scale_y(:);
%Dy = Dy0.*scale_y(:);

t_center = t(window_size/2:step:N-window_size/2);

reconstructed = zeros(size(t_center));
rec_no_offs = zeros(size(t_center));

w_best = w_best(:);
for i = 1:length(t_center)-1
    tt = t_center(i);     % Center time of window
    reconstructed(i) = Cy(i)*sin(w_best(i)*tt) + Dy(i)*cos(w_best(i)*tt) + By(i);
    rec_no_offs(i) = Cy(i)*sin(w_best(i)*tt) + Dy(i)*cos(w_best(i)*tt);
end

plot(t_center, reconstructed,'DisplayName','Reconstructed')
hold on 
plot(t_center, rec_no_offs,'DisplayName','Rec_no_offs')
plot(t,Ay, 'DisplayName','raw')
legend

%% phase

phi = atan2(Dy, Cy);

plot(t_center,phi)