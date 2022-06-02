function j_e = J_e(z, xi, xj)
j_e = zeros(3,6);
% dx/dxi(d는 편미분)
j_e(1,1) = -cos(xi(3)); j_e(1,2) = -sin(xi(3));
j_e(2,1) = sin(xi(3)); j_e(2,2) = -cos(xi(3));
j_e(1,3) = [-sin(xi(3)), cos(xi(3))]*(xj(1:2)-xi(1:2));
j_e(1,3) = [-cos(xi(3)), -sin(xi(3))]*(xj(1:2)-xi(1:2));
j_e(3,3) = -1;

% dx/dxj(d는 편미분)
j_e(1:2,4:5) = rot(xi(3));
j_e(3,6) = 1;
z_inv_tf = inv(v2t(z));
j_e = z_inv_tf*j_e;
end