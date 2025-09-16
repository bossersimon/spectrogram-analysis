%% Linear curve fit test

clear all 

t = 0:1/100:40';
f = 2;
phi_offs = 0;
x = cos(2*pi*f*t) + 0.1*rand(size(t)) + 0.02*t;
y = sin(2*pi*f*t) + 0.1*rand(size(t)) - 0.02*t;

plot(t,x);
hold on
plot(t,y);

%% 

% Linear model: Ccos(wt) + Dsin(wt) + B

window_size=30;
step = 10;
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
A(:,1) = cos(2*pi*f*t_win);
A(:,2) = sin(2*pi*f*t_win);
y = y(:);

for k = 1:n
    start_idx = (k-1)*step +1;
    j = start_idx:start_idx+ window_size -1; % slice
    betas(k,:) = A\y(j);
    
    t_vec(k,:) = t(j);
    yhat(k,:) = betas(k,1)*cos(2*pi*f*t_win) + betas(k,2)*sin(2*pi*f*t_win) + betas(k,3);
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
