

%% spectrogram

M = readmatrix("bil4.txt");

Mx = M(:,1);
My = M(:,2);

Mx2 = Mx;

Gz = M(:,6);


figure("Name","raw")

dt = 1/100;
N = size(M,1);
t = seconds((0:N-1)*dt);
tiledlayout(3,2,"TileIndexing","columnmajor")
ax = [];
for i = 1:6
    ax(end+1) =nexttile;
    plot(t, M(:,i))
    grid on
    % xticksformat("MM:SS")
end
linkaxes(ax,"x")


%%% DC removal
%for frame = 1:300:length(Mx)
%    Mx2(frame:frame+1) = Mx2(frame:frame+1);
%end

M2 = zeros(size(M));
for frame = 1:300:length(Mx)
    idx_end = min(frame + 300 - 1, length(Mx));
    current_window = M(frame:idx_end,1:3);
    current_window = current_window - mean(current_window);    
    M2(frame:idx_end,1:3) = current_window;
end
M2(:,4:6) = M(:,4:6);

%%
tiledlayout(3,2,"TileIndexing","columnmajor")
ax = [];
for i = 1:6
    ax(end+1) =nexttile;
    plot(t, M2(:,i))
    grid on
    % xticksformat("MM:SS")
end
linkaxes(ax,"x")

%%

fs = 100;
wsize = 300;
ovlap = 299;
Ndft = 1024;

w = rectwin(wsize); % window function can be changed to something else
[s,f,t] = spectrogram(Mx,w,ovlap,Ndft,fs);

figure;
spectrogram(Mx,w,ovlap,Ndft,fs);
view(90,-90)
%waterplot(s,f,t)

th = 0.5;   % threshold frequency in herz
th_idx = round(th/(fs/Ndft)) +1;   % index of roughly this frequency
threshold = (th_idx-1)*fs/Ndft;    % actual threshold with this index

[~,f0_relative_idx] = max(abs(s(th_idx:end,:)));  % returns relative index of masked array
f0_idx = f0_relative_idx + th_idx - 1;            % corresponding index in full array


f_vals = f(f0_idx)
f_avg = mean(f(f0_idx(30:end)))

wheel_radius = 2;
v_avg = f_avg*wheel_radius*3.6

cols = 1:size(s,2);   % cols for sub2ind
S = s(sub2ind(size(s), f0_idx, cols));  % values of S at wanted frequencies

tol = 1e-6;
S(abs(S) < tol) = 0;

phi_x = angle(S);
phix_unwrapped = unwrap(phi_x);
model_x = cos(f_avg*2*pi*t);
x_analytic = hilbert(model_x);
phi_est = unwrap(angle(x_analytic));
%phi_est = acos(model_x);

dt = 1/fs
gz = cumtrapz(Gz)*dt*pi/180;
t_gyro = (0:length(gz)-1) *dt;
gz_interp = interp1(t_gyro,gz,t,"linear","extrap");


figure;
plot(phix_unwrapped)
hold on
plot(-phix_unwrapped)
plot(model_x)
plot(phi_est)
plot(gz_interp)
legend('phi', '-phi','x(t)','phi_{model}','gz')
title("acc_x")
hold off
