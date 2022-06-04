function vec = t2v(tf)
vec = zeros(3,1);
vec(1) = tf(1,3);
vec(2) = tf(2,3);
vec(3) = atan2(tf(2,1), tf(1,1));
end