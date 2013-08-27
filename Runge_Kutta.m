function y = Runge_Kutta(f, t, h, y0)
% Runge Kutta 4th order method for solving ODE
% input: f = function handle
%        t = time span
%        h = time step

y = zeros(length(t), length(y0));
y(1,:) = y0;

for n = 1 : length(t) - 1
    k1 = f(t(n), y(n,:));
    k2 = f(t(n) + h/2, y(n,:) + h/2 * k1);
    k3 = f(t(n) + h/2, y(n,:) + h/2 * k2);
    k4 = f(t(n) + h,   y(n,:) + h   * k3);
    y(n+1,:) = y(n,:) + h/6*(k1+2*k2+2*k3+k4);
end