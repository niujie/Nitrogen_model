function [f] = f2(x)
% f2(x) = 0, with x = [x(1);x(2)]
% represents a system of 2 non-linear equations
f1 = x(1)^2 + x(2)^2 - 50;
f2 = x(1)*x(2) -25;
f = [f1;f2];
% end function