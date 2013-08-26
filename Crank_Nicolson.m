function Crank_Nicolson()

v = 10;
Diff = 0.5*v;

xl = 0;
xr = 50;
nx = 100;
tfinal = 8;

dx = (xr-xl)/nx;
dt = 0.08;
x = linspace(xl,xr,nx)';

nsteps = round(tfinal/dt);
nplot = 1;

tn = 0;
u0 = zeros(size(x));
u0(1) = 20.0;
u0(end) = 0;
u = u0;

ufine = u0(1)/2*(erfc((x-v*tn)/sqrt(4*Diff*tn))+exp(x*v/Diff).*erfc((x+v*tn)/sqrt(4*Diff*tn)));
plot(x,ufine,'r',x,u,'bo-')
drawnow

alpha = -v*dt/dx/4;
beta = Diff*dt/dx^2/2;
A = (alpha - beta) * ones(size(u));
A(1) = 0;
A(end) = 0;
B = (1 + 2*beta) * ones(size(u));
B(1) = 1;
B(end) = 1;
C = (-alpha - beta) * ones(size(u));
C(1) = 0;
C(end) = 0;
D = zeros(size(u));
I = (2:nx-1)';

for n = 1 : nsteps
    tnp = tn + dt;
    
    D(I) = (-alpha+beta)*u(I-1)+(1-2*beta)*u(I)+(alpha+beta)*u(I+1);
    D(1) = 20;
    D(end) = 0;
    
    u = TDMAsolver(A, B, C, D);
    if mod(n,nplot)==0 || n==nsteps
        ufine = u0(1)/2*(erfc((x-v*tnp)/sqrt(4*Diff*tnp))+exp(x*v/Diff).*erfc((x+v*tnp)/sqrt(4*Diff*tnp)));
        plot(x,ufine,'r',x,u,'bo-')
        title(sprintf('t = %9.5e  after %4i time steps with %5i grid points',...
                       tnp,n,nx))
        drawnow
    end
    
    tn = tnp;
end

end

%------------------------------------------------------------

function x = TDMAsolver(a,b,c,d)
%a, b, c are the column vectors for the compressed tridiagonal matrix, d is the right vector
n = length(d); % n is the number of rows
 
% Modify the first-row coefficients
c(1) = c(1) / b(1);    % Division by zero risk.
d(1) = d(1) / b(1);   
 
for i = 2:n-1
    temp = b(i) - a(i) * c(i-1);
    c(i) = c(i) / temp;
    d(i) = (d(i) - a(i) * d(i-1))/temp;
end
 
d(n) = (d(n) - a(n) * d(n-1))/( b(n) - a(n) * c(n-1));
 
% Now back substitute.
x(n) = d(n);
for i = n-1:-1:1
    x(i) = d(i) - c(i) * x(i + 1);
end
end