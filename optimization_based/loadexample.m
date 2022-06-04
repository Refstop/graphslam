function [z_1t, u_1t, c, landmark_n, T, ij] = loadexample()
x_dim = 3;
m_dim = 2;
var = 0.1;
landmark_n = 2;

u1 = [2,0,pi/2];
u2 = [2,0,pi/2];
u3 = [2,0,pi/2];
u4 = [2,0,0];
% u1 = [2,0,0];
% u2 = [0,2,0];
% u3 = [-2,0,0];
% u4 = [0,-2,0];
u_1t = [u1;u2;u3;u4]';

z1 = [-1,1;1,1];
z2 = [-1,1;-1,-1];
z3 = [-1,1;-3,1];
z4 = [-1,1;-1,3];
z_1t = {z1',z2',z3',z4'};

T = length(u_1t(1,:));
for i = 1:T
    ij(i,:) = [i,i+1];
end

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
end