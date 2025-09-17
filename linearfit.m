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

x = M_filt(:,1);
y = M_filt(:,2);

plot(t,x);
hold on
plot(t,y);
legend('x','y')

%% 

% Linear model: Ccos(wt) + Dsin(wt) + B

window_size=30;
step = 30;
N = length(y); % recording length

n = floor((N - window_size) / step) + 1; % number of fits

% preallocation
time = zeros(1,n);
%xvals = zeros(1,n);
%yvals = zeros(1,n);
betas = zeros(n,3);
yhat = zeros(n,window_size);
t_vec = zeros(n,window_size);

A = ones(window_size,3);     % design matrix
t_win = t(1:window_size);
y = y(:);

% Now need to iterate over different f

f0 = 0; % starting frequency
delta = 2; % 2*delta search range
best_y = zeros(n,window_size);

y_fit = zeros(1,window_size);
y_best = zeros(1,window_size);
    

for k = 1:n
    best_err = inf;
    best_f = 0;
    
    f_candidates = linspace(f0 - delta, f0 + delta, 300);

    for f = f_candidates
        A(:,1) = cos(2*pi*f*t_win);
        A(:,2) = sin(2*pi*f*t_win);

        start_idx = (k-1)*step +1;
        j = start_idx:start_idx+ window_size -1; % slice
        betas(k,:) = A\y(j);
        
        y_fit = betas(k,1)*cos(2*pi*f*t_win) + betas(k,2)*sin(2*pi*f*t_win) + betas(k,3);

        err = immse(y(j),y_fit);

        if err<best_err
            best_err = err;
            best_f = f;
            y_best = y_fit;
        end
    end

    f0 = best_f; % update search range

    yhat(k,:) = y_best; % store results
    t_vec(k,:) = t(j);

end 
    

%% plotting

figure;

plot(t,y);
hold on
for k=1:n
    start_idx = (k-1)*step +1;
    plt_indices = start_idx;
    plot(t_vec(k,:), yhat(k,:))
end

%plot(t,yhat,'r','LineWidth',2);

%[sin(w*t_rel), cos(w*t_rel), ones(size(t_rel))];

%t_windows = zeros(n,window_size);
%x_fits = zeros(n,window_size);
%y_fits = zeros(n,window_size);

%k=1;
%center_idx = window_size/2;
