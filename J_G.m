function G = J_G(mu,u,delta_t)
v_w = u(1)/u(2);
delta_th = u(2)*delta_t;
J_x = -v_w*cos(mu(3))+v_w*cos(mu(3)+delta_th);
J_y = -v_w*sin(mu(3))+v_w*sin(mu(3)+delta_th);
G = [1,0,J_x; 0,1,J_y; 0,0,1];
end