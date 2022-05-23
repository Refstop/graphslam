clc;clear all;close all;

%% GraphSLAM_initialize(u_1:t)
delta_t = 1;
mu = [0;0;0]; % x와 동일
u_1t = [1,realmin;1,realmin;1,realmin]; % 1s, 2s, 3s
for i = 1:length(u_1t)
    mu(:,i+1) = mu(:,i) + delta_xytheta(mu(:,i), u_1t(i,:), delta_t);
end
% plot(mu(1,:), mu(2,:), 'o')

%% GraphSLAM_linearize(u_1:t, z_1:t, c_1:t, mu_0:t)
Omega = 0; Xi = 0;
x0_info_mat = realmax*eye(3,3);
Omega = x0_info_mat;
R = 0.05*eye(3,3);
for i = 2:length(mu)
    x_hat = mu(:,i-1) + delta_xytheta(mu(:,i-1), u_1t(i-1,:), delta_t);
    G_t = J_G(mu(:,i-1), u_1t(i-1,:), delta_t);
    omega_temp = [-G_t';ones(3,3)]*inv(R)*[-G_t,ones(3,3)];
    Omega(pos3(i-1), pos3(i)) = omega_temp(1:3,1:3);
end