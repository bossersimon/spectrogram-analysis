
%% Test
clear all

t = (1:100)';
A = ones(100,3);
A(:,2) = cos((2*pi)/50*t);
A(:,3) = sin((2*pi)/50*t);
y = 2*cos((2*pi)/50*t+pi/2+randn(size(t)));
y = y(:);
beta = A\y;
yhat = beta(1)+beta(2)*cos((2*pi)/50*t) + beta(3)*sin((2*pi)/50*t);

plot(t,y,'b');
hold on
plot(t,yhat,'r','LineWidth',2);