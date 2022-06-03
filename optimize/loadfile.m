function [z_1t, u_1t, c, landmark_n, T, ij] = loadfile()
data = importdata("input.txt");
N = length(data.data);

measure = [];
t = [];
u_1t = [];
tmp = 0;
landmark_n = 1;
ij=[];
a=1;
for i = 1:N
    if strcmp(data.textdata{i}, 'EDGE_SE2')
        t = [t,data.data(i,2)];
        ij(a,:) = [a, a+1];
        a = a + 1;
        u_1t = [u_1t,data.data(i,3:5)'];
    elseif strcmp(data.textdata{i}, 'EDGE_SE2_XY')
        measure = [measure;data.data(i,1:4)];
    end
end
T = length(t);
c = cell(T,1);
z_1t = cell(T,1);
measure_N = length(measure);
for i = 1:measure_N
    if tmp ~= measure(i,1)
        tmp = measure(i,1);
        time = tmp - landmark_n + 1;
    end
    if tmp <= measure(i,2)
        fdict(landmark_n) = measure(i,2);
        landmark_n = landmark_n + 1;
    end
    c{time} = [c{time}, find(fdict == measure(i,2))];
    z_1t{time} = [z_1t{time}, measure(i,3:4)'];
end
landmark_n = landmark_n - 1;

end