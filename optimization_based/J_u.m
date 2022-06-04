function j_e = J_u(z, xi, xj)
j_e = zeros(3,6);
% dx/dxi(d는 편미분)
j_e(1,1) = -cos(xi(3)); j_e(1,2) = -sin(xi(3));
j_e(2,1) = sin(xi(3)); j_e(2,2) = -cos(xi(3));
j_e(1,3) = [-sin(xi(3)), cos(xi(3))]*(xj(1:2)-xi(1:2));
j_e(1,3) = [-cos(xi(3)), -sin(xi(3))]*(xj(1:2)-xi(1:2));
j_e(3,3) = -1;

% dx/dxj(d는 편미분)
j_e(1:2,4:5) = inv(rot(xi(3)));
j_e(3,6) = 1;
z_inv_tf = inv(v2t(z));
z_inv_tf(1:2,3) = 0;
j_e = z_inv_tf*j_e;

% rot_z = rot(z(3));
% ab_temp = rot_z*[-sin(xi(3))*cos(xj(3)) + cos(xi(3))*sin(xj(3));
%                 cos(xi(3))*cos(xj(3)) + sin(xi(3))*sin(xj(3))];
% a = ab_temp(1);
% b = ab_temp(2);
% x = a/b;
% 
% ab_p_temp = rot_z*[-cos(xi(3))*cos(xj(3)) - sin(xi(3))*sin(xj(3));
%                 -sin(xi(3))*cos(xj(3)) + cos(xi(3))*sin(xj(3))];
% a_p = ab_p_temp(1);
% b_p = ab_p_temp(2);
% j_e(3,3) = (a_p*b - a*b_p)/(b^2) * 1/(1+x^2);
% 
% ab_p_temp = rot_z*[sin(xi(3))*sin(xj(3)) + cos(xi(3))*cos(xj(3));
%                 -cos(xi(3))*sin(xj(3)) + sin(xi(3))*cos(xj(3))];
% a_p = ab_p_temp(1);
% b_p = ab_p_temp(2);
% j_e(3,6) = (a_p*b - a*b_p)/(b^2) * 1/(1+x^2);
% [j_e(3,3), j_e(3,6)]; % 왜 항상 -1, 1???
end