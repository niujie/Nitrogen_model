function [x,iter] = newtonm(x0,f,J) 
% Newton-Raphson method applied to a 
% system of linear equations f(x) = 0, 
% given the jacobian function J, with 
% J = del(f1,f2,...,fn)/del(x1,x2,...,xn) 
% x = [x1;x2;...;xn], f = [f1;f2;...;fn] 
% x0 is an initial guess of the solution 
N = 100; % define max. number of iterations 
epsilon = 1e-10; % define tolerance 
maxval = 10000.0; % define value for divergence 
xx = x0; % load initial guess 
dx = 0.001;
while (N>0) 
    JJ = feval(J,f,xx,dx); 
    if abs(det(JJ))<epsilon 
        error('newtonm - Jacobian is singular - try new x0'); 
    end; 
    xn = xx - JJ\feval(f,xx);
    if abs(feval(f,xn))<epsilon 
        x=xn; 
        iter = 100-N; 
        return; 
    end; 
    if abs(feval(f,xx))>maxval 
        iter = 100-N; 
        disp(['iterations = ',num2str(iter)]); 
        error('Solution diverges'); 
    end; 
    N = N - 1; 
    xx = xn; 
end; 
error('No convergence after 100 iterations.'); 
% end function