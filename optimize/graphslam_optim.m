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
for i = 1:T
    u_tf = u_tf*v2t(u_1t(:,i));
    x_plot(:,i+1) = t2v(u_tf);
end
plot(x_plot(1,:), x_plot(2,:),'o');
hold on;

%% test data from file
% [z_1t, u_1t, c, landmark_n, T, ij] = loadfile();
% mu = [0;0;0]; % x와 동일
% for t = 1:T
%     mu(:,t+1) = g(mu(:,t), u_1t(:,t));
% end
% plot(mu(1,:), mu(2,:), 'o')
% hold on;

%% graphslam using GN optimization
%% 1. pose graph
X = zeros(x_dim,T+1);

max_iter = 2;
for iter = 1:max_iter
    H = zeros(x_dim*(T+1));
    b = zeros(x_dim*(T+1),1);
    error_sum = 0;
    for i = 1:T
        id_i = ij(i,1); id_j = ij(i,2);
        x_i_tf = v2t(X(:,id_i));
        x_j_tf = v2t(X(:,id_j));
        z_ij_tf = v2t(u_1t(:,id_i));
        error = t2v(inv(z_ij_tf)*(inv(x_i_tf)*x_j_tf));
        Omega = R;
        J = J_e(u_1t(:,id_i), X(:,id_i), X(:,id_j));
        error_sum = error_sum + error'*error;
        H_ii = J(:,1:3)'*Omega*J(:,1:3);
        H_ij = J(:,1:3)'*Omega*J(:,4:6);
        H_jj = J(:,4:6)'*Omega*J(:,4:6);
        b_i = -J(:,1:3)'*Omega*error;
        b_j = -J(:,4:6)'*Omega*error;
        
        id_i_range = x_dim*(id_i-1)+1:x_dim*id_i;
        id_j_range = x_dim*(id_j-1)+1:x_dim*id_j;

        H(id_i_range, id_i_range) = H(id_i_range, id_i_range) + H_ii;
        H(id_j_range, id_j_range) = H(id_j_range, id_j_range) + H_jj;
        H(id_i_range, id_j_range) = H(id_i_range, id_j_range) + H_ij;
        H(id_j_range, id_i_range) = H(id_j_range, id_i_range) + H_ij';
        b(id_i_range) = b(id_i_range) + b_i;
        b(id_j_range) = b(id_j_range) + b_j;

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
    % ???
%     delta_x = inv(H)*b;
    X = X + reshape(delta_x,3,T+1);
end
plot(X(1,:), X(2,:));
xlim([-2 4]); ylim([-2 4]);
grid on;