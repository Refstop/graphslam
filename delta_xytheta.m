function delta_xyth = delta_xytheta(mu,u,delta_t)
v_w = u(1)/u(2);
delta_th = u(2)*delta_t;
delta_x = -v_w*sin(mu(3))+v_w*sin(mu(3)+delta_th);
delta_y = v_w*cos(mu(3))-v_w*cos(mu(3)+delta_th);
delta_xyth = [delta_x;delta_y;delta_th];
end