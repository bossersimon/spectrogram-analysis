

%% Mex spectogram stuff

M = readmatrix("bil1.txt");

Mx = M(:,1);
My = M(:,2);

fs = 100;
wsize = 500;
ovlap = 400;
Ndft = 1024;

w = rectwin(wsize);
%s = spectrogram(Mx);
%spectrogram(s(:,7),'yaxis')

[s,f,t] = spectrogram(Mx,w,ovlap,Ndft,fs);

waterplot(s,f,t)

th = 0.5;   % threshold frequency in herz
th_idx = round(th/(fs/Ndft)) +1;   % index of roughly this frequency
threshold = (th_idx-1)*fs/Ndft;    % actual threshold with this index

%phi = zeros(1,size(s,2)); 

[~,f0_relative_idx] = max(abs(s(th_idx:end,:)));  % returns relative index of masked array
f0_idx = f0_relative_idx + th_idx - 1;            % corresponding index in full array

cols = 1:size(s,2);   % cols for sub2ind
S = s(sub2ind(size(s), f0_idx, cols));  % values of S at wanted frequencies

phi_x = angle(S);

function waterplot(s,f,t)
% Waterfall plot of spectrogram
    waterfall(f,t,abs(s)'.^2)
    set(gca,XDir="reverse",View=[30 50])
    xlabel("Frequency (Hz)")
    ylabel("Time (s)")
end

