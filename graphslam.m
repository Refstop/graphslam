clc;clear all;close all;

%% GraphSLAM_initialize(u_1:t)
T = 3;
delta_t = 1;
landmark_n = 2;
mu = [0;0;0]; % x와 동일
for t = 1:T
    r = 0.1*randn(1,2);
    u_1t(t,:) = [1,0] + r; % 1s, 2s, 3s
end
z1 = [1,1,-pi/2;2,-1,pi/2];
z2 = [0,1,-pi/2;1,-1,pi/2];
z3 = [-1,1,-pi/2;0,-1,pi/2];
z4 = [-2,1,-pi/2;-1,-1,pi/2];
z_1t = [z1;z2;z3;z4];

for t = 1:T+1
    r = 0.1*randn(2,3);
    z_1t(pos(2,t), pos(3,1)) = z_1t(pos(2,t), pos(3,1)) + r;
end
for t = 1:T
    mu(:,t+1) = mu(:,t) + delta_xytheta(mu(:,t), u_1t(t,:), delta_t);
end
% plot(mu(1,:), mu(2,:), 'o')

%% GraphSLAM_linearize(u_1:t, z_1:t, c_1:t, mu_0:t)
Omega = realmax*eye(3,3); Xi = [0;0;0];
R = 0.05*eye(3,3);
for t = 2:length(mu)
    x_hat = mu(:,t-1) + delta_xytheta(mu(:,t-1), u_1t(t-1,:), delta_t);
    G_t = J_G(mu(:,t-1), u_1t(t-1,:), delta_t);
    omega_temp = [-G_t';ones(3,3)]*inv(R)*[-G_t,ones(3,3)];
    xi_temp = [-G_t';ones(3,3)]*inv(R)*[x_hat - G_t*mu(:,t-1)];
    Omega(pos(3,t), pos(3,t)) = zeros(3,3); Xi(pos(3,t)) = zeros(3,1);

    Omega(pos(3,t-1), pos(3,t-1)) = Omega(pos(3,t-1), pos(3,t-1)) + omega_temp(pos(3,1), pos(3,1));
    Omega(pos(3,t-1), pos(3,t)) = Omega(pos(3,t-1), pos(3,t)) + omega_temp(pos(3,1), pos(3,2));
    Omega(pos(3,t), pos(3,t-1)) = Omega(pos(3,t), pos(3,t-1)) + omega_temp(pos(3,2), pos(3,1));
    Omega(pos(3,t), pos(3,t)) = Omega(pos(3,t), pos(3,t)) + omega_temp(pos(3,2), pos(3,2));
    
    Xi(pos(3,t-1)) = Xi(pos(3,t-1)) + xi_temp(pos(3,1));
    Xi(pos(3,t)) = Xi(pos(3,t)) + xi_temp(pos(3,2));
end
mu_t = mu(:,T+1);
z_1t = xy2rb(z_1t, mu);
map = zeros(3*landmark_n, 1);
for t = 1:1 % T+1
    Q = 0.05*eye(3,3);
    z = z_1t(pos(2,t),:);
    for j = 1:landmark_n
        map_k = map(pos(3,j));
        delta = map_k(1:2) - mu(1:2,t);
        sqrt_q = norm(delta);
        z_hat = [sqrt_q; atan2(delta(1),delta(2))-mu(3,t); map_k(3)-mu(3,t)];
        H = (1/(sqrt_q^2))*[-sqrt_q*delta(1), -sqrt_q*delta(2), 0, sqrt_q*delta(1), sqrt_q*delta(2), 0;
            delta(2), -delta(1), -sqrt_q^2, -delta(2), delta(1), 0;
            0, 0, -sqrt_q^2, 0, 0, sqrt_q^2];
        Omega(pos(3,T+1+j), pos(3,T+1+j)) = zeros(3,3); Xi(pos(3,T+1+j)) = zeros(3,1);
        omega_temp = H'*inv(Q)*H; xi_temp = H'*inv(Q)*(z(j,:)' - z_hat + H*[mu(:,t);map_k]);

        Omega(pos(3,t), pos(3,t)) = Omega(pos(3,t), pos(3,t)) + omega_temp(pos(3,1), pos(3,1));
        Omega(pos(3,t), pos(3,T+1+j)) = Omega(pos(3,t), pos(3,T+1+j)) + omega_temp(pos(3,1), pos(3,2));
        Omega(pos(3,T+1+j), pos(3,t)) = Omega(pos(3,T+1+j), pos(3,t)) + omega_temp(pos(3,2), pos(3,1));
        Omega(pos(3,T+1+j), pos(3,T+1+j)) = Omega(pos(3,T+1+j), pos(3,T+1+j)) + omega_temp(pos(3,2), pos(3,2));

        Xi(pos(3,t)) = Xi(pos(3,t)) + xi_temp(pos(3,1));
        Xi(pos(3,T+1+j)) = Xi(pos(3,T+1+j)) + xi_temp(pos(3,2));
    end
end