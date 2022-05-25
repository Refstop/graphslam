clc;clear all;close all;

%% GraphSLAM_initialize(u_1:t)
T = 3;
delta_t = 1;
landmark_n = 2;
mu = [0;0;0]; % x와 동일
for i = 1:T
    r = 0.1*randn(1,2);
    u_1t(i,:) = [1,0] + r; % 1s, 2s, 3s
end
z1 = [1,1,-pi/2;2,-1,pi/2];
z2 = [0,1,-pi/2;1,-1,pi/2];
z3 = [-1,1,-pi/2;0,-1,pi/2];
z4 = [-2,1,-pi/2;-1,-1,pi/2];
z_1t = [z1;z2;z3;z4];

for i = 1:T+1
    r = 0.1*randn(2,3);
    z_1t(pos(2,i), pos(3,1)) = z_1t(pos(2,i), pos(3,1)) + r;
end
for i = 1:T
    mu(:,i+1) = mu(:,i) + delta_xytheta(mu(:,i), u_1t(i,:), delta_t);
end
% plot(mu(1,:), mu(2,:), 'o')

%% GraphSLAM_linearize(u_1:t, z_1:t, c_1:t, mu_0:t)
Omega = realmax*eye(3,3); Xi = [0;0;0];
R = 0.05*eye(3,3);
for i = 2:length(mu)
    x_hat = mu(:,i-1) + delta_xytheta(mu(:,i-1), u_1t(i-1,:), delta_t);
    G_t = J_G(mu(:,i-1), u_1t(i-1,:), delta_t);
    omega_temp = [-G_t';ones(3,3)]*inv(R)*[-G_t,ones(3,3)];
    xi_temp = [-G_t';ones(3,3)]*inv(R)*[x_hat - G_t*mu(:,i-1)];
    Omega(pos(3,i), pos(3,i)) = zeros(3,3); Xi(pos(3,i)) = zeros(3,1);
%     for a = 1:2
%         Xi(pos(3,i-2+a)) = Xi(pos(3,i-2+a)) + xi_temp(pos(3,a));
%         for b = 1:2
%             Omega(pos(3,i-2+a), pos(3,i-2+b)) = Omega(pos(3,i-2+a), pos(3,i-2+b)) + omega_temp(pos(3,a), pos(3,b));
%         end
%     end
    Omega(pos(3,i-1), pos(3,i-1)) = Omega(pos(3,i-1), pos(3,i-1)) + omega_temp(pos(3,1), pos(3,1));
    Omega(pos(3,i-1), pos(3,i)) = Omega(pos(3,i-1), pos(3,i)) + omega_temp(pos(3,1), pos(3,2));
    Omega(pos(3,i), pos(3,i-1)) = Omega(pos(3,i), pos(3,i-1)) + omega_temp(pos(3,2), pos(3,1));
    Omega(pos(3,i), pos(3,i)) = Omega(pos(3,i), pos(3,i)) + omega_temp(pos(3,2), pos(3,2));
    
    Xi(pos(3,i-1)) = Xi(pos(3,i-1)) + xi_temp(pos(3,1));
    Xi(pos(3,i)) = Xi(pos(3,i)) + xi_temp(pos(3,2));
end
mu_t = mu(:,T+1);
z_1t = xy2rb(z_1t, mu);
map = zeros(3*landmark_n, 1);
for i = 1:T
    Q = 0.05*eye(3,3);
    for j = 1:T+1
        for k = 1:landmark_n
            map_k = map(pos(3,k));
            delta = map_k(1:2) - mu(1:2,i);
            q = norm(delta);
            z_hat = [q; atan2(delta(1),delta(2))-mu(3,i); map_k(3)-mu(3,i)];
        end
    end
end