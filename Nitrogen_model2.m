function Nitrogen_model2()
% parameters and model equations refer to
% Mark A. Widdowson et al., A Numerical Transport Model for Oxygen-
% and Nitrate-Based Respiration Linked to Substrate and Nutritent
% Availability in Porous Media. WRR 1988

v       = 10;        % unit: cm/day, average linear velocity
Diff    = 0.5*v;     % unit: cm2/day, longitudinal dispersivity (cm) * velocity (cm/day)
xl      = 0;         % unit: cm, left boundary
xr      = 50;        % unit: cm, right boundary
nx      = 100;       % number of nodes
tfinal  = 8;         % unit: day, total simulation time
dx      = (xr-xl)/nx;% unit: cm, space step
dt      = 0.008;      % unit: day, time step
x       = linspace(xl, xr, nx)';  % unit: cm, space grids
nsteps  = round(tfinal/dt);       % time steps
nplot   = 1;

tn      = 0;         % t(n)
ns      = 9;         % number of different component
u0      = zeros(nx, ns);
% initial values, unit: mg/cm3
% S = u0(:,1)   % substrate concentration in pore fluid
% O = u0(:,2)   % oxygen    concentration in pore fluid
% N = u0(:,3)   % nitrate   concentratoin in pore fluid
% A = u0(:,4)   % ammonia   concentration in pore fluid
% M = u0(:,5)   % biomass   concentration
% s = u0(:,6)   % substrate concentration within colony
% o = u0(:,7)   % oxygen    concentration within colony
% n = u0(:,8)   % nitrate   concentration within colony
% a = u0(:,9)   % ammonia   concentration within colony
u0(:,1) = 1.0e-3;
u0(:,2) = 2.0e-3;
u0(:,3) = 0.0e-3;
u0(:,4) = 1.0e-3;
u0(:,5) = 5.65e-4;
u0(:,6) = u0(:,1)/10;
u0(:,7) = u0(:,2)/10;
u0(:,8) = u0(:,3)/10;
u0(:,9) = 1.0e-3;
% left boundary values, unit: mg/cm3
u0(1,1) = 10.0e-3;
u0(1,2) = 2.0e-3;
u0(1,4) = 1.0e-3;
u0(1,5) = 5.65e-4;
% right boundary values
u0(nx,1) = 0.0;
u0(nx,2) = 0.0;
u0(nx,3) = 0.0;
u0(nx,4) = 1.0e-3;
u0(nx,5) = 0.0;
u0(nx,6) = 0.0;
u0(nx,7) = 0.0;
u0(nx,8) = 0.0;
u0(nx,9) = 1.0e-3;

u = u0;

% construct left hand matrix
alpha  = -v*dt/dx/4;
beta   = Diff*dt/dx^2/2;
A      = (alpha - beta) * ones(nx);
A(1)   = 0;
A(nx) = 0;
B      = (1 + 2*beta) * ones(nx);
B(1)   = 1;
B(nx) = 1;
C      = (-alpha - beta) * ones(nx);
C(1)   = 0;
C(nx) = 0;

I = (2:nx-1)';
% legend text for plot
specs  = {'SUB', 'O2', 'NO3', 'NH4', 'Bio', 'sub', 'o2', 'no3', 'nh4'};
% plot initial condition in the unit of mg/l
plot(x, u*1e3);
legend(specs(1:ns));
title(sprintf('t = %6.2f  after %4i time steps with %5i grid points',...
               tn,0,nx))
xlim([0 xr])
xlabel('DISTANCE (cm)')
ylabel('CONCENTRATION (mg/l)')
drawnow    

for n = 1 : nsteps
    tnp = tn + dt;
    
    % transport for S, O, N, A
    for i = 1 : 4
        D      = zeros(nx);
        D(I)   = (-alpha+beta)*u(I-1,i)+(1-2*beta)*u(I,i)+(alpha+beta)*u(I+1,i);
        D(1)   = u0(1,i);
        D(nx)  = 2*u(nx-1,i) - u(nx-2,i);
        u(:,i) = TDMAsolver(A, B, C, D);
    end
    u(1,5:9) = u0(1,5:9);
    u(nx,5:9) = u0(nx,5:9);
    for i = 1 : nx
        error = 1e8;
        uold = u(i,:);
        while error > 1e-5
            unew = u(i,:);
            % Newton iterative method for s, o, n, a ==> u(:,6:9)
            % u(:,1:4) are S, O, N, A        
            u(i,6:9) = newtonm(uold(6:9)',u(i,1:4)',@fun,@jacobFD);
            % Runge Kutta method solving ODE equation (25) for biomass, u(:,5)
            % note that M is update based on the old step value and the
            % iterated values of s, o, n and a, not the itrated value of M
            %u(i,5) = Runge_Kutta_mode(@mode, tnp, dt, uold(5), u(i,6:9));
            u(i,5) = mode_tp(dt, uold(5), u(i,6:9));
            % Runge Kutta method solving system of ODEs in eq (9) to (12) for
            % S, O, N, and A
            for item = 1 : 4
                % note S, O, N and A are updated based on the ADE solution
                % not on the previous iterated solution
                u(i,item) = Runge_Kutta_odes(@odes, tnp, dt, uold(item), u(i,item+5), u(i,5), item);
                %u(i,item) = odes_tp(dt, uold(item), u(i,item+5), u(i,5), item);
            end
            error = max(abs(u(i,:)-unew));
        end
    end
    
    if mod(n,nplot)==0 || n==nsteps
        flag = 1:ns;
        for i = 1 : ns
            if all(u(:,i) < 1e-8)
                flag(i) = 0;
            end
        end
        flag(flag==0) = [];
        plot(x, u(:,flag)*1e3)
        legend(specs(flag))
        title(sprintf('t = %6.2f  after %4i time steps with %5i grid points',...
                       tnp,n,nx))
        xlim([0 50])
        ylim([0 10])
        xlabel('DISTANCE (cm)')
        ylabel('CONCENTRATION (mg/l)')
        drawnow
    end    
    tn = tnp;    
end

end


%--------------------------------------------------------------------------

function y = Runge_Kutta_mode(f, t, h, y, u)
% Runge Kutta 4th order method for solving ODE
% input: f = function handle
%        h = time step

k1 = f(t, y, u);
k2 = f(t + h/2, y' + h/2 * k1, u);
k3 = f(t + h/2, y' + h/2 * k2, u);
k4 = f(t + h,   y' + h   * k3, u);
y = y + h/6*(k1'+2*k2'+2*k3'+k4');

end

%--------------------------------------------------------------------------
function dMdt = mode(~, M, u)
% equatio (25) in the paper
% all parameters should be the same to those in fun.m

muo     = 3.1;        % unit: day-1, maximum specific growth rate for the active heterotrophic population
mun     = 2.9;        % unit: day-1, maximum specific growth rate
ko      = 0.02;       % unit: day-1, microbial decay coefficient for aerobic respiration
kn      = 0.02;       % unit: day-1, microbial decay coefficient for nitrate respiration
Kso     = 0.04;       % unit: mg/cm3, substrate saturation constants under aerobic respiration
Ksn     = 0.04;       % unit: mg/cm3, substrate saturation constants under nitrate respiration
Ko      = 0.00077;    % unit: mg/cm3, oxygen saturation constants under aerobic respiration
Kn      = 0.0026;     % unit: mg/cm3, nitrate saturation constants under nitrate respiration
Kao     = 0.001;      % unit: mg/cm3, ammonia saturation constants under aerobic respiration
Kan     = 0.001;      % unit: mg/cm3, ammonia saturation constants under nitrate respiration
Kc      = 1;          % unit: mass/length3, inhibition coefficient

s = u(1);
o = u(2);
n = u(3);
a = u(4);

%dMdt = M * ((muo*s*o*a/(Kso+s)/(Ko+o)/(Kao+a)-ko) + ...
%    (mun*s*n*a/(Ksn+s)/(Kn+n)/(Kan+a)-kn)/(1+o/Kc));
% ammonia is in excess
dMdt = M * ((muo*s*o/(Kso+s)/(Ko+o)-ko) + ...
    (mun*s*n/(Ksn+s)/(Kn+n)-kn)/(1+o/Kc));

end


%--------------------------------------------------------------------------
function M = mode_tp(dt, M, u)
% equatio (25) in the paper
% all parameters should be the same to those in fun.m
% trapezoid method to implicitly solve ODE

muo     = 3.1;        % unit: day-1, maximum specific growth rate for the active heterotrophic population
mun     = 2.9;        % unit: day-1, maximum specific growth rate
ko      = 0.02;       % unit: day-1, microbial decay coefficient for aerobic respiration
kn      = 0.02;       % unit: day-1, microbial decay coefficient for nitrate respiration
Kso     = 0.04;       % unit: mg/cm3, substrate saturation constants under aerobic respiration
Ksn     = 0.04;       % unit: mg/cm3, substrate saturation constants under nitrate respiration
Ko      = 0.00077;    % unit: mg/cm3, oxygen saturation constants under aerobic respiration
Kn      = 0.0026;     % unit: mg/cm3, nitrate saturation constants under nitrate respiration
Kao     = 0.001;      % unit: mg/cm3, ammonia saturation constants under aerobic respiration
Kan     = 0.001;      % unit: mg/cm3, ammonia saturation constants under nitrate respiration
Kc      = 1;          % unit: mass/length3, inhibition coefficient

s = u(1);
o = u(2);
n = u(3);
a = u(4);

%dMdt = M * ((muo*s*o*a/(Kso+s)/(Ko+o)/(Kao+a)-ko) + ...
%    (mun*s*n*a/(Ksn+s)/(Kn+n)/(Kan+a)-kn)/(1+o/Kc));
% ammonia is in excess
f = (muo*s*o/(Kso+s)/(Ko+o)-ko) + ...
    (mun*s*n/(Ksn+s)/(Kn+n)-kn)/(1+o/Kc);

M = (1+dt*f/2)/(1-dt*f/2)*M;

end
%--------------------------------------------------------------------------

function y = Runge_Kutta_odes(f, t, h, y, u, M, item)
% Runge Kutta 4th order method for solving ODE
% input: f = function handle
%        h = time step

k1 = f(t, y, u, M, item);
k2 = f(t + h/2, y' + h/2 * k1, u, M, item);
k3 = f(t + h/2, y' + h/2 * k2, u, M, item);
k4 = f(t + h,   y' + h   * k3, u, M, item);
y = y + h/6*(k1'+2*k2'+2*k3'+k4');

end

%--------------------------------------------------------------------------
function dUdt = odes(~, U, u, M, item)
% ODE system in equations (9) to (12) in the paper
% inputs:
%         U is one of [S, O, N, A],   u(:,1:4)
%         u is one of [s, o, n, a],   u(:,6:9)
%         M is biomass concentration, u(:,5)
%         item is an integer to specify which one in u(:,1:4)

% all parameters should be the same to those in fun.m
Dsb     = 1.03;       % unit: cm2/day, coefficient of diffusion in the boundary layer
Dob     = 2.19;       % unit: cm2/day, coefficient of diffusion in the boundary layer
Dnb     = 1.50;       % unit: cm2/day, coefficient of diffusion in the boundary layer
Dab     = 1.86;       % unit: cm2/day, coefficient of diffusion in the boundary layer
mc      = 2.83e-8;    % unit: mg/colony, cell mass of an average colony
delta   = 50.0e-4;    % unit: um -> cm, boundary layer thickness
rc      = 10.0e-4;    % unit: um -> cm, disk of uniform radius attached to aquifer sediment
tau     = 1.0e-4;     % unit: um -> cm, disk thickness attached to aquifer sediment
beta    = pi*rc^2 + 2*pi*rc*tau;    % total area across which diffusion occurs

fs      = 1.10;       % substrate retardatoin factor (ratio of actual pore 
                      % water velocity to velocity of the species in solution)
fo      = 1.00;       % oxygen retardatoin factor
fn      = 1.00;       % nitrate retardatoin factor
fa      = 2.20;       % ammonia retardatoin factor
theta   = 1.0;        % porosity

switch item
    case 1
        D = Dsb;      % for substrate
        f = fs;
    case 2
        D = Dob;      % for Oxygen
        f = fo;
    case 3
        D = Dnb;      % for nitrate
        f = fn;
    case 4
        D = Dab;      % for ammonia
        f = fa;
    otherwise
        error('number of item exceeds the maximum 4')
end

dUdt = -D*M/theta*(U-u)/delta*beta/mc/f;

end

%--------------------------------------------------------------------------
function U = odes_tp(dt, U, u, M, item)
% ODE system in equations (9) to (12) in the paper
% trapezoid method to implicitly solve ODEs
% inputs:
%         U is one of [S, O, N, A],   u(:,1:4)
%         u is one of [s, o, n, a],   u(:,6:9)
%         M is biomass concentration, u(:,5)
%         item is an integer to specify which one in u(:,1:4)

% all parameters should be the same to those in fun.m
Dsb     = 1.03;       % unit: cm2/day, coefficient of diffusion in the boundary layer
Dob     = 2.19;       % unit: cm2/day, coefficient of diffusion in the boundary layer
Dnb     = 1.50;       % unit: cm2/day, coefficient of diffusion in the boundary layer
Dab     = 1.86;       % unit: cm2/day, coefficient of diffusion in the boundary layer
mc      = 2.83e-8;    % unit: mg/colony, cell mass of an average colony
delta   = 50.0e-4;    % unit: um -> cm, boundary layer thickness
rc      = 10.0e-4;    % unit: um -> cm, disk of uniform radius attached to aquifer sediment
tau     = 1.0e-4;     % unit: um -> cm, disk thickness attached to aquifer sediment
beta    = pi*rc^2 + 2*pi*rc*tau;    % total area across which diffusion occurs

fs      = 1.10;       % substrate retardatoin factor (ratio of actual pore 
                      % water velocity to velocity of the species in solution)
fo      = 1.00;       % oxygen retardatoin factor
fn      = 1.00;       % nitrate retardatoin factor
fa      = 2.20;       % ammonia retardatoin factor
theta   = 1.0;        % porosity

switch item
    case 1
        D = Dsb;      % for substrate
        f = fs;
    case 2
        D = Dob;      % for Oxygen
        f = fo;
    case 3
        D = Dnb;      % for nitrate
        f = fn;
    case 4
        D = Dab;      % for ammonia
        f = fa;
    otherwise
        error('number of item exceeds the maximum 4')
end

coef = -D*M/theta/delta*beta/mc/f;

U = ((1+dt*coef/2)*U - dt*coef*u)/(1-dt*coef/2);

end