function j_e = J_z(x, m)
j_e = zeros(2,5);
j_e(:,1:2) = -rot(x(3));
P = [-sin(x(3)), -cos(x(3));
    cos(x(3)), -sin(x(3))];
j_e(:,3) = P * (m(1:2)-x(1:2));
j_e(:,4:5) = rot(x(3));
end