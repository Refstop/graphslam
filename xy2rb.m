function z_rb = xy2rb(z)
r = sqrt(z(:,1).^2 + z(:,2).^2);

phi = atan2(z(:,2), z(:,1));
s = z(:,3);
z_rb = [r, phi, s];
end