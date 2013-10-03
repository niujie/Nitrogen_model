function Nitrogen_model()
% parameters and model equations refer to
% Mee-Sun Lee et al., Nitrogen transformation and transport modeling
% in groundwater aquifers. Ecological Modelling, 2006

v       = 10;        % unit: cm/day, average linear velocity
Diff    = 0.5*v;     % unit: cm2/day, longitudinal dispersivity (cm) * velocity (cm/day)
xl      = 0;         % unit: cm, left boundary
xr      = 50;        % unit: cm, right boundary
nx      = 100;       % number of nodes
tfinal  = 4;         % unit: day, total simulation time
dx      = (xr-xl)/nx;% unit: cm, space step
dt      = 0.08;      % unit: day, time step
x       = linspace(xl, xr, nx)';  % unit: cm, space grids
nsteps  = round(tfinal/dt);       % time steps
nplot   = 1;

tn      = 0;         % t(n)
ns      = 9;         % number of different component
u0      = zeros(nx, ns);
% initial values, unit: mg/l
% NH4  = u0(:,1)
% NO2  = u0(:,2)
% NO3  = u0(:,3)
% N2   = u0(:,4)
% CH2O = u0(:,5)
% O2   = u0(:,6)
% X1   = u0(:,7)
% X2   = u0(:,8)
% X3   = u0(:,9)
u0(:,1) = 1.0;
u0(:,2) = 0;
u0(:,3) = 5.0;
u0(:,4) = 0.0;
u0(:,5) = 1.0;
u0(:,6) = 2.0;
u0(:,7) = 0;
u0(:,8) = 0;
u0(:,9) = 3.64;
% left boundary values, unit: mg/l
u0(1,1) = 1.0;
u0(1,2) = 0;
u0(1,3) = 5.0;
u0(1,4) = 0.0;
u0(1,5) = 20.0;
u0(1,6) = 2.0;
u0(1,7) = 0;
u0(1,8) = 0;
u0(1,9) = 3.64;
% right boundary values
u0(nx,1) = 0.0;
u0(nx,2) = 0.0;
u0(nx,3) = 0.0;
u0(nx,4) = 0.0;
u0(nx,5) = 0.0;
u0(nx,6) = 0.0;
u0(nx,7) = 0.0;
u0(nx,8) = 0.0;
u0(nx,9) = 0.0;

u = u0;

I = (2:nx-1)';
% legend text for plot
specs  = {'NH4', 'NO2', 'NO3', 'N2', 'CH2O', 'O2', 'X1', 'X2', 'X3'};
% plot initial condition
plot(x, u);
legend(specs(1:ns));
title(sprintf('t = %6.2f  after %4i time steps with %5i grid points',...
               tn,0,nx))
xlim([0 xr])
xlabel('DISTANCE (cm)')
ylabel('CONCENTRATION (mg/l)')
drawnow

for n = 1 : nsteps    
    tnp = tn + dt;
    
    % transport for each component
    for i = 1 : ns
        % construct left hand matrix
        ne           = 0.44;    % dimensionless, porosity
        rho          = 1.56e6;  % unit: mg/l, soil bulk density
        KdCH2O       = 2.93e-8; % unit: l/mg, linear partitioning coefficient of DOC
        KdNH4        = 35.2e-8; % unit: l/mg, linear partitioning coefficient of ammonium
        switch i
            case {2, 3, 4, 6, 7, 8, 9}
                R = 1;
            case 1
                R = 1 + rho/ne * KdNH4;
            case 5
                R = 1 + rho/ne * KdCH2O;
            otherwise
        end
        alpha  = -v*dt/dx/4/R;
        beta   = Diff*dt/dx^2/2/R;
        A      = (alpha - beta) * ones(nx);
        A(1)   = 0;
        A(nx) = 0;
        B      = (1 + 2*beta) * ones(nx);
        B(1)   = 1;
        B(nx) = 1;
        C      = (-alpha - beta) * ones(nx);
        C(1)   = 0;
        C(nx) = 0;        
        D      = zeros(nx);
        D(I)   = (-alpha+beta)*u(I-1,i)+(1-2*beta)*u(I,i)+(alpha+beta)*u(I+1,i);
        D(1)   = u0(1,i);
        D(nx)  = 2*u(nx-1,i) - u(nx-2,i);
        u(:,i) = TDMAsolver(A, B, C, D);
    end
    
    % Runge Kutta method solving system of ODEs
    for i = 1 : nx
        u(i,:) = Runge_Kutta(@myode, tnp, dt, u(i,:));
    end
    
    if mod(n,nplot)==0 || n==nsteps
        flag = 1:ns;
        for i = 1 : ns
            if all(u(:,i) < 1e-3)
                flag(i) = 0;
            end
        end
        flag(flag==0) = [];
        plot(x, u(:,flag))
        legend(specs(flag))
        title(sprintf('t = %6.2f  after %4i time steps with %5i grid points',...
                       tnp,n,nx))
        xlim([0 50])
        xlabel('DISTANCE (cm)')
        ylabel('CONCENTRATION (mg/l)')
        drawnow
    end    
    tn = tnp;    
end

end


%--------------------------------------------------------------------------

function y = Runge_Kutta(f, t, h, y)
% Runge Kutta 4th order method for solving ODE
% input: f = function handle
%        h = time step

k1 = f(t, y);
k2 = f(t + h/2, y + h/2 * k1);
k3 = f(t + h/2, y + h/2 * k2);
k4 = f(t + h,   y + h   * k3);
y = y + h/6*(k1+2*k2+2*k3+k4);

end

%--------------------------------------------------------------------------
function dydt = myode(~,y)

mu_max_nit1  = 10;      % unit: 1/day, the maximum substrate utilization rate for nitrification (ammonium -> nitrite)
mu_max_nit2  = 10;      % unit: 1/day, the maximum substrate utilization rate for nitrification (nitrite -> nitrate)
mu_max_denit = 40;      % unit: 1/day, the maximum substrate utilization rate for denitrification (nitrate -> nitrogen)
mu_max_oxid  = 30;      % unit: 1/day, the maximum substrate utilization rate for organic carbon oxidation
K_NH4        = 1;       % unit: mg/l, half saturation constant for NH4+
K_NO2        = 1.8;     % unit: mg/l, half saturation constant for NO2-
K_NO3        = 2.6;     % unit: mg/l, half saturation constant for NO3-
K_O2         = 0.77;    % unit: mg/l, half saturation constant for O2
K_CH2O       = 40;      % unit: mg/l, half saturation constant for CH2O
k_O2I        = 0.01;    % unit: mg/l, oxygen inhibition constant
kb1          = 1;       % unit: mg/l, ammonia-oxidizing biomass inhibition constant
kb2          = 1;       % unit: mg/l, nitrite-oxidizing biomass inhibition constant
kb3          = 0.5;     % unit: mg/l, heterotrophic biomass inhibition constant
yNO2_NH4     = 2.5504;  % ratio of NH4  substrate to NO2 substrate consumed
yNO3_NO2     = 1.3478;  % ratio of NO2  substrate to NO3 substrate consumed
yN2_NO3      = 0.4518;  % ratio of NO3  substrate to N2  substrate consumed
yO2_NH4      = 1.7739;  % ratio of NH4  substrate to O2  substrate consumed
yO2_NO2      = 0.6955;  % ratio of NO2  substrate to O2  substrate consumed
yO2_CH2O     = 1.0657;  % ratio of CH2O substrate to O2  substrate consumed
Y1           = 0.45;    % microbial yield coefficients for autotrophic ammonia-oxidizing bacteria
Y2           = 0.45;    % microbial yield coefficients for autotrophic nitrite-oxidizing bacteria
Y3           = 0.5;     % microbial yield coefficients for heterotrophic bacteria
d1           = 0.02;    % unit: 1/day, death or maintenance rate constant of autotrophic ammonia-oxidizing bacteria
d2           = 0.02;    % unit: 1/day, death or maintenance rate constant of autotrophic nitrite-oxidizing bacteria
d3           = 0.02;    % unit: 1/day, death or maintenance rate constant of heterotrophic biomass
ne           = 0.44;    % dimensionless, porosity
rho          = 1.56e6;  % unit: mg/l, soil bulk density
KdCH2O       = 2.93e-8; % unit: l/mg, linear partitioning coefficient of DOC
KdNH4        = 35.2e-8; % unit: l/mg, linear partitioning coefficient of ammonium
R_NH4        = 1 + rho/ne * KdNH4;       % retardation factor for NH4+
R_NO2        = 1;       % retardation factor for NO2-
R_NO3        = 1;       % retardation factor for NO3-
R_N2         = 1;       % retardation factor for N2
R_CH2O       = 1 + rho/ne * KdCH2O;       % retardation factor for CH2O
R_O2         = 1;       % retardation factor for O2

NH4  = y(1);
NO2  = y(2);
NO3  = y(3);
%N2   = y(4);
CH2O = y(5);
O2   = y(6);
X1   = y(7);    % unit: mg/l, concentration of autotrophic ammonia-oxidizing biomass
X2   = y(8);    % unit: mg/l, concentration of autotrophic nitrite-oxidizing biomass
X3   = y(9);    % unit: mg/l, the heterotrophic biomass concentration

r1 = mu_max_nit1*X1*(kb1/(kb1+X1))*(NH4/(K_NH4+NH4))*(O2/(K_O2+O2));    % unit: mg/day, substrate utilization rate by ammonium oxidation process 1
%r1 = mu_max_nit1*X1*(kb1/(kb1+X1))*1.0*(O2/(K_O2+O2));                  % ammonia is in excess
r2 = mu_max_nit2*X2*(kb2/(kb2+X2))*(NO2/(K_NO2+NO2))*(O2/(K_O2+O2));    % unit: mg/day, substrate utilization rate by nitrite oxidation process 2
r3 = mu_max_denit*X3*(kb3/(kb3+X3))*(k_O2I/(k_O2I+O2))*(CH2O/(K_CH2O+CH2O))*(NO3/(K_NO3+NO3));  % unit: mg/day, substrate utilization rate by denitrification
r4 = mu_max_oxid*X3*(kb3/(kb3+X3))*(CH2O/(K_CH2O+CH2O))*(O2/(K_O2+O2)); % unit: mg/day, substrate utilization rate by denitrification

dNH4dt  = - r1 / R_NH4;
dNO2dt  = (yNO2_NH4 * r1 - r2) / R_NO2;
dNO3dt  = (yNO3_NO2 * r2 - r3) / R_NO3;
dN2dt   = (yN2_NO3 * r3) / R_N2;
dCH2Odt = - r4 / R_CH2O;
dO2dt   = (-yO2_NH4 * r1 - yO2_NO2 * r2 - yO2_CH2O * r4) / R_O2;
dX1dt   = Y1 * r1 - X1 * d1;
dX2dt   = Y2 * r2 - X2 * d2;
dX3dt   = Y3 * (r3+r4) - X3 * d3;

dydt = zeros(1,9);
dydt(1) = dNH4dt;
dydt(2) = dNO2dt;
dydt(3) = dNO3dt;
dydt(4) = dN2dt;
dydt(5) = dCH2Odt;
dydt(6) = dO2dt;
dydt(7) = dX1dt;
dydt(8) = dX2dt;
dydt(9) = dX3dt;

end