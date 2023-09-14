%% simple random walk
close; clear; clc;

p = 0.5;
N = 10000;
Xi = 0;
X = [];
X = [X, Xi];

for i = 1 : N
    Zi = binornd(1,p);
    if Zi == 0
        Zi = -1;
    end
    Xi = Xi + Zi;
    X = [X, Xi];
end

plot(X);