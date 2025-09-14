
%% Test
clear all

t = (1:100)';
X = ones(100,3);
X(:,2) = cos((2*pi)/50*t);
X(:,3) = sin((2*pi)/50*t);
y = 2*cos((2*pi)/50*t-pi/4+randn(size(t)));
y = y(:);
beta = X\y;
yhat = beta(1)+beta(2)*cos((2*pi)/50*t) + beta(3)*sin((2*pi)/50*t);

plot(t,y,'b');
hold on
plot(t,yhat,'r','LineWidth',2);