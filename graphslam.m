clc;clear all;close all;
%% params
x_dim = 3;
m_dim = 2;
var = 0.05;
mu = [0;0;0]; % x와 동일

%% test example

% landmark_n = 2;
% 
% u1 = [2,0,0];
% u2 = [0,2,pi/2];
% u3 = [0,2,pi/2];
% u4 = [0,2,pi/2];
% u_1t = [u1;u2;u3;u4]';
% 
% z1 = [-1,1;1,1];
% z2 = [-1,1;-1,-1];
% z3 = [-1,1;-3,1];
% z4 = [-1,1;-1,3];
% z_1t = {z1',z2',z3',z4'};
% 
% T = length(u_1t(1,:));
% 
% for t = 1:T
%     r = var*randn(x_dim,1);
%     u_1t(:,t) = u_1t(:,t) + r; % 1s, 2s, 3s
% end
% 
% 
% for t = 1:T
%     r = var*randn(landmark_n,m_dim);
%     z_1t{t} = z_1t{t} + r;
% end
% 
% c1 = [1,2];
% c2 = [1,2];
% c3 = [1,2];
% c4 = [1,2];
% c = {c1;c2;c3;c4};

%% load data from file
data = importdata("input.txt");
N = length(data.data);

measure = [];
t = [];
u_1t = [];
tmp = 0;
landmark_n = 1;

for i = 1:N
    if strcmp(data.textdata{i}, 'EDGE_SE2')
        t = [t,data.data(i,2)];
        u_1t = [u_1t,data.data(i,3:5)'];
    elseif strcmp(data.textdata{i}, 'EDGE_SE2_XY')
        measure = [measure;data.data(i,1:4)];
    end
end
T = length(t);
c = cell(T,1);
z_1t = cell(T,1);
measure_N = length(measure);
for i = 1:measure_N
    if tmp ~= measure(i,1)
        tmp = measure(i,1);
        time = tmp - landmark_n + 1;
    end
    if tmp <= measure(i,2)
        fdict(landmark_n) = measure(i,2);
        landmark_n = landmark_n + 1;
    end
    c{time} = [c{time}, find(fdict == measure(i,2))];
    z_1t{time} = [z_1t{time}, measure(i,3:4)'];
end
landmark_n = landmark_n - 1;

%% GraphSLAM_initialize(u_1:t)
for t = 1:T
    mu(:,t+1) = g(mu(:,t), u_1t(:,t));
end

%% GraphSLAM_linearize(u_1:t, z_1:t, c_1:t, mu_0:t)
size = x_dim*(T+1) + m_dim*landmark_n;
Omega = zeros(size, size); Xi = zeros(size, 1);
Omega(1:3,1:3) = realmax*eye(3,3); Xi(1:3) = [0;0;0];
R = 0.1*eye(3,3);
R(3,3) = deg2rad(1);
for t = 2:length(mu)
    x_hat = g(mu(:,t-1), u_1t(:,t-1));
    G_t = J_G(mu(:,t-1), u_1t(:,t-1));
    A = [-G_t, eye(3,3)];
    omega_temp = [-G_t';eye(3,3)]*inv(R)*[-G_t,eye(3,3)];
    xi_temp = [-G_t';eye(3,3)]*inv(R)*[x_hat - G_t*mu(:,t-1)];
    
    Omega(x_dim*(t-2)+1:x_dim*(t-2)+6, x_dim*(t-2)+1:x_dim*(t-2)+6) = Omega(x_dim*(t-2)+1:x_dim*(t-2)+6, x_dim*(t-2)+1:x_dim*(t-2)+6) + omega_temp;
    Xi(x_dim*(t-2)+1:x_dim*(t-2)+6) = Xi(x_dim*(t-2)+1:x_dim*(t-2)+6) + xi_temp;
end

for t = 1:T % T+1
    Q = 0.1*eye(m_dim,m_dim);
    z = z_1t{t};
    c_N = length(c{t});
    for i = 1:c_N
        j = c{t}(i);
        mu_j = mu(1:2,t+1) + rot(mu(3,t+1))*z(:,i);
        z_hat = h(mu(:,t+1),mu_j);
        H = J_H(mu(:,t+1),mu_j);

        pos_x = x_dim*t;
        pos_m = x_dim*(T+1) + m_dim*(j-1);

        omega_temp = H'*inv(Q)*H;
        xi_temp = H'*inv(Q)*(z(:,i) - z_hat + H*[mu(:,t+1);mu_j]);

        Omega(pos_x+1:pos_x+x_dim, pos_x+1:pos_x+x_dim) = Omega(pos_x+1:pos_x+x_dim, pos_x+1:pos_x+x_dim) + omega_temp(1:3, 1:3);
        Omega(pos_x+1:pos_x+x_dim, pos_m+1:pos_m+m_dim) = Omega(pos_x+1:pos_x+x_dim, pos_m+1:pos_m+m_dim) + omega_temp(1:3, 4:5);
        Omega(pos_m+1:pos_m+m_dim, pos_x+1:pos_x+x_dim) = Omega(pos_m+1:pos_m+m_dim, pos_x+1:pos_x+x_dim) + omega_temp(4:5, 1:3);
        Omega(pos_m+1:pos_m+m_dim, pos_m+1:pos_m+m_dim) = Omega(pos_m+1:pos_m+m_dim, pos_m+1:pos_m+m_dim) + omega_temp(4:5, 4:5);

        Xi(pos_x+1:pos_x+x_dim) = Xi(pos_x+1:pos_x+x_dim) + xi_temp(1:3);
        Xi(pos_m+1:pos_m+m_dim) = Xi(pos_m+1:pos_m+m_dim) + xi_temp(4:5);
    end
end

%% GraphSLAM_reduce(Omega, Xi)
Omega_til = Omega(1:3*(T+1),1:3*(T+1));
Xi_til = Xi(1:3*(T+1));
Omega_til = Omega_til - Omega(1:3*(T+1),3*(T+1)+1:end)*inv(Omega(3*(T+1)+1:end,3*(T+1)+1:end))*Omega(3*(T+1)+1:end,1:3*(T+1));
Xi_til = Xi_til - Omega(1:3*(T+1),3*(T+1)+1:end)*inv(Omega(3*(T+1)+1:end,3*(T+1)+1:end))*Xi(3*(T+1)+1:end);

%% GraphSLAM_solve(Omega_til, Xi_til)
om_dim = x_dim*(T+1);
mu_final = zeros(x_dim*(T+1) + m_dim*landmark_n, 1);
Sigma = inv(Omega_til);
mu_final(1:om_dim) = Sigma*Xi_til;
Sigma_m = inv(Omega(om_dim+1:end, om_dim+1:end));
mu_final(om_dim+1:end) = Sigma_m * (Xi(om_dim+1:end) - Omega(om_dim+1:end, 1:om_dim)*mu_final(1:om_dim));

%% plot
for i = 1:T+1
    mu_f_x(i) = mu_final(3*(i-1)+1);
    mu_f_y(i) = mu_final(3*(i-1)+2);
end
for i = 1:landmark_n
    mu_f_lm_x(i) = mu_final(3*(T+1)+2*(i-1)+1);
    mu_f_lm_y(i) = mu_final(3*(T+1)+2*(i-1)+2);
end
plot(mu(1,:), mu(2,:), 'o')
% xlim([-1 4]); ylim([-1.2 1.2]);
hold on;
plot(mu_f_x,mu_f_y)
hold on;
plot(mu_f_lm_x,mu_f_lm_y,'o')
grid on;