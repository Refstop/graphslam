clc;clear all;close all;
%% params
x_dim = 3;
m_dim = 2;
var = 0.1;
mu = [0;0;0]; % x와 동일

R = var*eye(x_dim,x_dim);
R(3,3) = deg2rad(1);

Q = var*eye(m_dim,m_dim);

%% test data
[z_1t, u_1t, c, landmark_n, T, ij] = loadexample();
x_plot = [0;0;0];
u_tf = eye(3);
for t = 1:T
    u_tf = u_tf*v2t(u_1t(:,t));
    x_plot(:,t+1) = t2v(u_tf);
end
plot(x_plot(1,:), x_plot(2,:),'o');
grid on;
hold on;

%% test data from file
% [z_1t, u_1t, c, landmark_n, T, ij] = loadfile();
% x_plot = [0;0;0];
% u_tf = eye(3);
% for t = 1:T
%     u_tf = u_tf*v2t(u_1t(:,t));
%     x_plot(:,t+1) = t2v(u_tf);
% end
% plot(x_plot(1,:), x_plot(2,:),'o');
% grid on;
% hold on;

%% graphslam using GN optimization
X = zeros(x_dim + m_dim*landmark_n, T+1);
Omega_u = inv(R);
Omega_z = inv(Q);

max_iter = 10;
for iter = 1:max_iter
    H = zeros(x_dim*(T+1) + m_dim*landmark_n*T);
    b = zeros(x_dim*(T+1) + m_dim*landmark_n*T,1);
    error_sum = 0;
    for t = 1:T
        % pose graph H,b
        id_i = ij(t,1); id_j = ij(t,2);
        x_i_tf = v2t(X(1:3,id_i));
        x_j_tf = v2t(X(1:3,id_j));
        z_ij_tf = v2t(u_1t(:,id_i));
        error_u = t2v(inv(z_ij_tf)*(inv(x_i_tf)*x_j_tf));
        
        J_ue = J_u(u_1t(:,id_i), X(1:3,id_i), X(1:3,id_j));
        error_sum = error_sum + error_u'*error_u;
        H_ii = J_ue(:,1:3)'*Omega_u*J_ue(:,1:3);
        H_jj = J_ue(:,4:6)'*Omega_u*J_ue(:,4:6);
        H_ij = J_ue(:,1:3)'*Omega_u*J_ue(:,4:6);
        b_i = -J_ue(:,1:3)'*Omega_u*error_u;
        b_j = -J_ue(:,4:6)'*Omega_u*error_u;
        
        id_i_range = x_dim*(id_i-1)+1:x_dim*id_i;
        id_j_range = x_dim*(id_j-1)+1:x_dim*id_j;

        H(id_i_range, id_i_range) = H(id_i_range, id_i_range) + H_ii;
        H(id_j_range, id_j_range) = H(id_j_range, id_j_range) + H_jj;
        H(id_i_range, id_j_range) = H(id_i_range, id_j_range) + H_ij;
        H(id_j_range, id_i_range) = H(id_j_range, id_i_range) + H_ij';
        b(id_i_range) = b(id_i_range) + b_i;
        b(id_j_range) = b(id_j_range) + b_j;

        % landmark information 추가
        z = z_1t{t};
        c_N = length(c{t});
        for i = 1:c_N
            j = c{t}(i);
            m_range = x_dim + m_dim*(i-1)+1:x_dim + m_dim*i;
            error_z = h(X(1:3,t+1), X(m_range,t)) - z(:,j);
            J_ze = J_z(X(1:3,t+1), X(m_range,t));
            H_ii = J_ze(:,1:3)'*Omega_z*J_ze(:,1:3);
            H_jj = J_ze(:,4:5)'*Omega_z*J_ze(:,4:5);
            H_ij = J_ze(:,1:3)'*Omega_z*J_ze(:,4:5);
            b_i = -J_ze(:,1:3)'*Omega_z*error_z;
            b_j = -J_ze(:,4:5)'*Omega_z*error_z;

            id_x_range = x_dim*(t-1)+1:x_dim*t;
            id_m_range = x_dim*(T+1) + m_dim*landmark_n*(t-1) + m_dim*(j-1)+1:x_dim*(T+1) + m_dim*landmark_n*(t-1) + m_dim*j;

            H(id_x_range, id_x_range) = H(id_x_range, id_x_range) + H_ii;
            H(id_m_range, id_m_range) = H(id_m_range, id_m_range) + H_jj;
            H(id_x_range, id_m_range) = H(id_x_range, id_m_range) + H_ij;
            H(id_m_range, id_x_range) = H(id_m_range, id_x_range) + H_ij';
            b(id_x_range) = b(id_x_range) + b_i;
            b(id_m_range) = b(id_m_range) + b_j;
        end
        % 바로바로 대입이 더 빠름
%         F_J = zeros(x_dim, x_dim*(T+1));
%         F_J(:,id_i_range) = J(:,1:3);
%         F_J(:,id_j_range) = J(:,4:6);
%         H = H + F_J'*Omega*F_J;
%         b = b + -F_J'*Omega*error;
    end
    % ???
    H(1:3,1:3) = H(1:3,1:3) + eye(3); % lm?
    SH=sparse(H);
    delta_x=SH\b;
    % delta_x = inv(H)*b;
    X(1:x_dim,:) = X(1:x_dim,:) + reshape(delta_x(1:x_dim*(T+1)), x_dim, T+1);
    X(x_dim+1:end,2:end) = X(x_dim+1:end,2:end) + reshape(delta_x(x_dim*(T+1)+1:end), m_dim*landmark_n, T);
end
plot(X(1,:), X(2,:));
% xlim([-2 4]); ylim([-2 4]);
grid on;
hold on;
plot(X(4,T+1), X(5,T+1), 'o');
hold on;
plot(X(6,T+1), X(7,T+1), 'o');
hold on;