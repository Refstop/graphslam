function z_rb = xy2rb(z, mu)
r = sqrt(z(:,1).^2 + z(:,2).^2);
mu_temp = mu(3,:);
for i = 1:length(mu_temp)
    mu_th(pos(2,i), 1) = [mu_temp(i); mu_temp(i)];
end
phi = atan2(z(:,2), z(:,1)) - mu_th;
s = z(:,3);
z_rb = [r, phi, s];
end