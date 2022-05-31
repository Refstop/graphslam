clc;clear all;close all;

%% create test dataset
x = [-5:2:10]';
% x = x + randn(1,length(x));
y = x.^2;
plot(x,y,'o');
hold on;
grid on;

%% Gauss-newton
p = ones(3,1);

max_iter = 10;


for i = 1:max_iter
    error = y - (p(1)*x.^2+p(2)*x+p(3));
    j_e = 2*p(1)*x+p(2);
    p = p - inv(j_e'*j_e)*j_e'*error;
end
