%% Linear curve fit test

clear all 


t = 0:1/100:40;
f = ones(1,length(t));
f = 2*f;  % 2Hz
phi_offs = 0;
Ax = cos(2*pi*f.*t) + 0.1*rand(size(t)) + 0.02*t;
Ay = sin(2*pi*f.*t) + 0.1*rand(size(t)) - 0.02*t;

plot(t,Ax);
hold on
plot(t,Ay);


%% 


% Linear model: Ccos(wt) + Dsin(wt) + B

window_size=30;
step = 10;
N = length(Ay); % recording length

n = floor((N - window_size) / step) + 1; % number of fits

% preallocation
n = floor((N - window_size) / step) + 1;
time = zeros(1,n);
xvals = zeros(1,n);
yvals = zeros(1,n);

A = ones(n,3);     % design matrix  
A(:,1) = cos(2*pi*f*t);
A(:,2) = sin(2*pi*f*t);

%[sin(w*t_rel), cos(w*t_rel), ones(size(t_rel))];

%t_windows = zeros(n,window_size);
%x_fits = zeros(n,window_size);
%y_fits = zeros(n,window_size);

%k=1;
%center_idx = window_size/2;

for w = 1:n
    
end
