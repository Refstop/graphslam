clc;clear all;close all;

%% create test dataset
x = [-5:2:10]';
% x = x + randn(1,length(x));
y = 2*x.^2;
plot(x,y,'o');
hold on;
grid on;

%% Gauss-newton
% p(1)*x.^2+p(2)*x+p(3) 모델이 필요 - 모델이 왜 필요해?
p = ones(3,1);

max_iter = 10;

for i = 1:max_iter
    error = y - (p(1)*x.^2+p(2)*x+p(3)); % 에러값이 필요
    j_e(:,1) = -x.^2;
    j_e(:,2) = -x;
    j_e(:,3) = -1; % 자코비안이 필요
    p = p - inv(j_e'*j_e)*j_e'*error;
end
x_a = [-5:0.1:10]';
y_a = p(1)*x_a.^2+p(2)*x_a+p(3);
plot(x_a,y_a);