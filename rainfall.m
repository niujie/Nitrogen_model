function p = rainfall()

p = 0.0;
n   = 0.4;          % Dimensionless, porosity
Zr  = 0.8;          % unit: m, soil depth
times = poissrnd(0.23);
for i = 1 : times
    p = p + exprnd(11/(n*Zr))*(n*Zr);
end
% pd1 = makedist('Poisson', 'lambda', 0.23);
% pd2 = makedist('Exponential', 'mu', 11);
% times = random(pd1);
% for i = 1 : times
%     p = p + random(pd2);
% end
p = p / 1000;   % mm to m