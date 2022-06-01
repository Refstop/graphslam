function j_e = J_e(z, xi, xj)
j_e = zeros(3,6);
% dx/dxi(d는 편미분)
j_e(1,1) = -cos(z(3))*cos(xi(3)) - sin(z(3))*sin(xi(3));
j_e(1,2) = cos(z(3))*sin(xi(3)) - sin(z(3))*cos(xi(3));
j_e(2,1) = sin(z(3))*cos(xi(3)) - cos(z(3))*sin(xi(3));
j_e(2,2) = -sin(z(3))*sin(xi(3)) - cos(z(3))*cos(xi(3));
a = -sin(xi(3))*xj(1) + cos(xi(3))*xj(2) + sin(xi(3))*xi(1) + cos(xi(3))*xi(2);
b = -cos(xi(3))*xj(1) - sin(xi(3))*xj(2) - cos(xi(3))*xi(1) + sin(xi(3))*xi(2);
j_e(1,3) = cos(z(3))*a + sin(z(3))*b;
j_e(2,3) = cos(z(3))*a + sin(z(3))*b;
rot_z = rot(z(3));
ab_temp = rot_z*[-sin(xi(3))*cos(xj(3)) + cos(xi(3))*sin(xj(3));
                cos(xi(3))*cos(xj(3)) + sin(xi(3))*sin(xj(3))];
a = ab_temp(1);
b = ab_temp(2);
x = a/b;

ab_p_temp = rot_z*[-cos(xi(3))*cos(xj(3)) - sin(xi(3))*sin(xj(3));
                -sin(xi(3))*cos(xj(3)) + cos(xi(3))*sin(xj(3))];
a_p = ab_p_temp(1);
b_p = ab_p_temp(2);
j_e(3,3) = (a_p*b - a*b_p)/(b^2) * 1/(1+x^2);

% dx/dxj(d는 편미분)
j_e(1,4) = cos(z(3))*cos(xi(3)) - sin(z(3))*sin(xi(3));
j_e(1,5) = cos(z(3))*sin(xi(3)) + sin(z(3))*cos(xi(3));
j_e(2,4) = -sin(z(3))*cos(xi(3)) - cos(z(3))*sin(xi(3));
j_e(2,5) = -sin(z(3))*sin(xi(3)) + cos(z(3))*cos(xi(3));

ab_p_temp = rot_z*[sin(xi(3))*sin(xj(3)) + cos(xi(3))*cos(xj(3));
                -cos(xi(3))*sin(xj(3)) + sin(xi(3))*cos(xj(3))];
a_p = ab_p_temp(1);
b_p = ab_p_temp(2);
j_e(3,6) = (a_p*b - a*b_p)/(b^2) * 1/(1+x^2);
end