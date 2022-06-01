function tf = v2t(vec)
tf = [cos(vec(3)), -sin(vec(3)), vec(1);
    sin(vec(3)), cos(vec(3)), vec(2);
    0, 0, 1];
end