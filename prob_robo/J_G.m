function G = J_G(x,u)
J_x = -u(1)*sin(x(3))-u(2)*cos(x(3));
J_y = u(1)*cos(x(3))-u(2)*sin(x(3));
G = [1,0,J_x; 0,1,J_y; 0,0,1];
end