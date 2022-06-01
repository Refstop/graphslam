clc;clear all;close all;
%% params
x_dim = 3;
m_dim = 2;
var = 0.05;
mu = [0;0;0]; % x와 동일

R = 0.1*eye(x_dim,x_dim);
R(3,3) = deg2rad(1);

Q = 0.1*eye(m_dim,m_dim);

%% test code
landmark_n = 2;

u1 = [2,0,0];
u2 = [0,2,pi/2];
u3 = [0,2,pi/2];
u4 = [0,2,pi/2];
u_1t = [u1;u2;u3;u4]';

z1 = [-1,1;1,1];
z2 = [-1,1;-1,-1];
z3 = [-1,1;-3,1];
z4 = [-1,1;-1,3];
z_1t = {z1',z2',z3',z4'};

T = length(u_1t(1,:));

for t = 1:T
    r = var*randn(x_dim,1);
    u_1t(:,t) = u_1t(:,t) + r; % 1s, 2s, 3s
end


for t = 1:T
    r = var*randn(landmark_n,m_dim);
    z_1t{t} = z_1t{t} + r;
end

c1 = [1,2];
c2 = [1,2];
c3 = [1,2];
c4 = [1,2];
c = {c1;c2;c3;c4};
%% test plot
X_plot = [0;0;0];
for i = 1:T
    X_plot(:,i+1) = v2t(u_1t(:,i))*X_plot(:,i)
end
plot(X_plot(1,:))

%% graphslam using GN optimization
%% 1. pose graph - z값 안씀
X = zeros(x_dim,T+1);
error_sum = 0;
H = 0; b = 0;
max_iter = 100;
for iter = 1:max_iter
    for i = 1:T
        x_i_tf = v2t(X(:,i));
        x_j_tf = v2t(X(:,i+1));
        z_ij_tf = v2t(u_1t(:,i));
        error = t2v(inv(z_ij_tf)*(inv(x_i_tf)*x_j_tf));
        Omega = R;
        error_sum = error_sum + error'*Omega*error;
        J = J_e(u_1t(:,i), X(:,i+1), X(:,i));
        F_J = zeros(x_dim, x_dim*(T+1));
        F_J(:,x_dim*(i-1)+1:x_dim*(i-1)+3) = J(:,1:3);
        F_J(:,x_dim*i+1:x_dim*i+3) = J(:,4:6);
        b = b + error'*Omega*F_J;
        H = H + F_J'*Omega*F_J;
    end
    H(H<1e-10) = 0;
    delta_x = -inv(H)*b';
    for i = 1:T+1
        X(:,i) = X(:,i) + delta_x(x_dim*(i-1)+1:x_dim*(i-1)+3);
    end
end
plot(X(1,:), X(2,:));
grid on;