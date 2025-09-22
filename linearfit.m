%% Linear curve fit test

clear all 
% 
% t = 0:1/100:40';
% f1 = 2;
% phi_offs = 0;
%x = cos(2*pi*f1*t) + 0.1*rand(size(t)) + 0.02*t;
%y = sin(2*pi*f1*t) + 0.1*rand(size(t)) - 0.02*t;

%x = chirp(t,0,40,5)+0.05*rand(size(t)) + 0.02*t;
%y = chirp(t,0,40,5)+0.05*rand(size(t)) - 0.02*t;

accelScale = 1/9.82; % scale accelerometer readings


M = readmatrix("recordings/recording_20250701_02.csv");

% Lowpass 
fc = 6; 
fs = 100;
n = 100; % filter order
b = fir1(n, (fc/(fs/2)), 'low');

M_filt = M;
M_filt(:,1:3) = filtfilt(b,1,M(:,1:3));

dt = 1/100;
N = size(M,1);
t = (0:N-1)*dt;
t= transpose(t);

%x = M_filt(:,1);
%y = M_filt(:,2);

x = M(:,1);
y = M(:,2);

plot(t,x);
hold on
plot(t,y);
legend('x','y')

%% 

% Linear model: Ccos(wt) + Dsin(wt) + B

window_size=80;
step = 80;
N = length(y); % recording length

n = floor((N - window_size) / step) + 1; % number of fits

% preallocation
time = zeros(1,n);
betas_y = zeros(n,3);
betas_x = zeros(n,3);
yhat = zeros(n,window_size);
xhat = zeros(n,window_size);
t_vec = zeros(n,window_size);

A = ones(window_size,3);     % design matrix
t_win = t(1:window_size);
y = y(:);

% Now need to iterate over different f

f0 = 0; % starting frequency
delta = 1; % 2*delta search range
best_y = zeros(n,window_size);
best_x = zeros(n,window_size);

y_fit = zeros(1,window_size);
y_best = zeros(1,window_size);
    
x_fit = zeros(1,window_size);
x_best = zeros(1,window_size);

tic
for k = 1:n
    best_err_y = inf;
    best_f_y = 0;

    best_err_x = inf;
    best_f_x = 0;
    
    f_candidates = linspace(f0 - delta, f0 + delta, 20);

    for f = f_candidates
        A(:,1) = cos(2*pi*f*t_win);
        A(:,2) = sin(2*pi*f*t_win);

        start_idx = (k-1)*step +1;
        j = start_idx:start_idx+ window_size -1; % slice
        betas_y(k,:) = A\y(j);
        betas_x(k,:) = A\x(j);
        
        y_fit = betas_y(k,1)*cos(2*pi*f*t_win) + betas_y(k,2)*sin(2*pi*f*t_win) + betas_y(k,3);
        x_fit = betas_x(k,1)*cos(2*pi*f*t_win) + betas_x(k,2)*sin(2*pi*f*t_win) + betas_x(k,3);

        err_y = immse(y(j),y_fit);
        err_x = immse(x(j),x_fit);

        if err_y<best_err_y
            best_err_y = err_y;
            best_f_y = f;
            y_best = y_fit;% -betas_y(k,3);
        end
        
        if err_x<best_err_x
            best_err_x = err_x;
            best_f_x = f;
            x_best = x_fit;% -betas_x(k,3);
        end
    end

    f0 = (best_f_y+best_f_x)/2; % update search range

    yhat(k,:) = y_best; % store results
    xhat(k,:) = x_best; % store results
    t_vec(k,:) = t(j);

end 
toc

%% plotting


fig = figure('Units','normalized','OuterPosition',[0 0 1 1]); 
set(fig, 'PaperOrientation', 'landscape');

plot(t,y,'Color','k');
hold on
%plot(t,x);
for k=1:n
    start_idx = (k-1)*step +1;
    plt_indices = start_idx;
    plot(t_vec(k,:), yhat(k,:)-betas_y(k,3), 'LineWidth',2)
    %plot(t_vec(k,:), xhat(k,:))
end
grid on

%%
fig = gcf;
exportgraphics(fig, 'linearfit_01.pdf', 'ContentType', 'vector');

%% phase

args = zeros(n,window_size);

yhat0 = yhat-betas_y(:,3);
xhat0 = xhat-betas_x(:,3);

figure;
for k=1:n
    start_idx = (k-1)*step +1;
    plt_indices = start_idx;

    args(k,:) = atan2(yhat0(k,:),xhat0(k,:));
    plot(t_vec(k,:), args(k,:), 'LineWidth',2)
    hold on
end

%plot(t,yhat,'r','LineWidth',2);

%[sin(w*t_rel), cos(w*t_rel), ones(size(t_rel))];

%t_windows = zeros(n,window_size);
%x_fits = zeros(n,window_size);
%y_fits = zeros(n,window_size);

%k=1;
%center_idx = window_size/2;
