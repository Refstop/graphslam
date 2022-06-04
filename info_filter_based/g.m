function g = g(x,u)
g = zeros(3,1);
g(1:2) = x(1:2) + rot(x(3))*u(1:2);
g(3) = x(3) + u(3);
end