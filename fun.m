function [f, varargout] = fun(x, vals)
% nonlinear system of equations (20) - (23)
% from Mark A. Widdowson et al., A Numerical Transport Model for Oxygen-
% and Nitrate-Based Respiration Linked to Substrate and Nutritent
% Availability in Porous Media. WRR 1988
% Outpus: function values and Jacobian matrix of the function
% Inputs
s       = x(1);       % unit: mass/length3, concentration of substrate within colony
o       = x(2);       % unit: mass/length3, concentration of oxygen    within colony
n       = x(3);       % unit: mass/length3, concentration of nitrate   within colony
a       = x(4);       % unit: mass/length3, concentration of ammonia   within colony
S       = vals(1);    % unit: mass/length3, pore fluid concentration of substrate
O       = vals(2);    % unit: mass/length3, pore fluid concentration of oxygen
N       = vals(3);    % unit: mass/length3, pore fluid concentration of nitrate
A       = vals(4);    % unit: mass/length3, pore fluid concentration of ammonia

% parameters
%rho     = 90.0;       % unit: mg/cm3, cell mass density
mc      = 2.83e-8;    % unit: mg/colony, cell mass of an average colony
%M       = 5.65e-4;    % unit: mg/cm3, biomass concentration
delta   = 50.0e-4;    % unit: um -> cm, boundary layer thickness
rc      = 10.0e-4;    % unit: um -> cm, disk of uniform radius attached to aquifer sediment
tau     = 1.0e-4;     % unit: um -> cm, disk thickness attached to aquifer sediment
beta    = pi*rc^2 + 2*pi*rc*tau;    % total area across which diffusion occurs

Dsb     = 1.03;       % unit: cm2/day, coefficient of diffusion in the boundary layer
Dob     = 2.19;       % unit: cm2/day, coefficient of diffusion in the boundary layer
Dnb     = 1.50;       % unit: cm2/day, coefficient of diffusion in the boundary layer
Dab     = 1.86;       % unit: cm2/day, coefficient of diffusion in the boundary layer
muo     = 3.1;        % unit: day-1, maximum specific growth rate for the active heterotrophic population
mun     = 2.9;        % unit: day-1, maximum specific growth rate
Yo      = 0.45;       % heterotrophic yield coefficient
Yn      = 0.5;        % heterotrophic yield coefficient
alphao  = 0.0402;     % the oxygen use coefficient for energy of maintenance for heterotrophic microorganisms
alphan  = 0.1;        % nitrate use coefficient for enery of maintenance for heterotropic microorganisms
gamma   = 1.4;        % oxygen use coefficient for synthesis of heterotrophic biomass
eta     = 2.2;        % nitrate use coefficient for synthesis of heterotrophic biomass
psi     = 0.122;      % ammonia use coefficient for productoin of heterotrophic biomass under aerobic respiration
epsilon = 0.122;      % ammonia use coefficient for productoin of heterotrophic biomass under nitrate respiration
ko      = 0.02;       % unit: day-1, microbial decay coefficient for aerobic respiration
kn      = 0.02;       % unit: day-1, microbial decay coefficient for nitrate respiration
Kso     = 0.04;       % unit: mg/cm3, substrate saturation constants under aerobic respiration
Ksn     = 0.04;       % unit: mg/cm3, substrate saturation constants under nitrate respiration
Ko      = 0.00077;    % unit: mg/cm3, oxygen saturation constants under aerobic respiration
Kn      = 0.0026;     % unit: mg/cm3, nitrate saturation constants under nitrate respiration
Kao     = 0.001;      % unit: mg/cm3, ammonia saturation constants under aerobic respiration
Kan     = 0.001;      % unit: mg/cm3, ammonia saturation constants under nitrate respiration
Kop     = 0.00077;    % unit: mg/cm3, oxygen saturation constant for decay
Knp     = 0.0026;     % unit: mg/cm3, nitrate saturaton constant for decay
Kc      = 1;          % unit: mass/length3, inhibition coefficient


f = zeros(size(x));   % (1) == eq(20), (2) = eq(21), (3) = eq(22), (4) = eq(23)
f(1) = -Dsb*(S-s)*beta/delta/mc + muo/Yo*s*o*a/(Kso+s)/(Ko+o)/(Kao+a) + ...
    mun/Yn*s*n*a/(Ksn+s)/(Kn+n)/(Kan+a)/(1+o/Kc);
f(2) = -Dob*(O-o)*beta/delta/mc + gamma*muo*s*o*a/(Kso+s)/(Ko+o)/(Kao+a) + ...
    alphao*ko*o/(Kop+o);
f(3) = -Dnb*(N-n)*beta/delta/mc + (eta*mun*s*n*a/(Ksn+s)/(Kn+n)/(Kan+a) + ...
    alphan*kn*n/(Knp+n))/(1+o/Kc);
f(4) = -Dab*(A-a)*beta/delta/mc + psi*muo*s*o*a/(Kso+s)/(Ko+o)/(Kao+a) + ...
    epsilon*mun*s*n*a/(Ksn+s)/(Kn+n)/(Kan+a)/(1+o/Kc);

% calculate f prime
if nargout > 1
    fp = zeros(size(x,1), size(x,1));

    fp(1,1) = (Dsb*beta)/(delta*mc) + (a*muo*o)/(Yo*(Kao + a)*(Ko + o)*(Kso + s)) ...
        - (a*muo*o*s)/(Yo*(Kao + a)*(Ko + o)*(Kso + s)^2) + ...
        (a*mun*n)/(Yn*(Kan + a)*(Kn + n)*(Ksn + s)*(o/Kc + 1)) - ...
        (a*mun*n*s)/(Yn*(Kan + a)*(Kn + n)*(Ksn + s)^2*(o/Kc + 1));
    fp(1,2) = (a*muo*s)/(Yo*(Kao + a)*(Ko + o)*(Kso + s)) - ...
        (a*muo*o*s)/(Yo*(Kao + a)*(Ko + o)^2*(Kso + s)) - ...
        (a*mun*n*s)/(Kc*Yn*(Kan + a)*(Kn + n)*(Ksn + s)*(o/Kc + 1)^2);
    fp(1,3) = (a*mun*s)/(Yn*(Kan + a)*(Kn + n)*(Ksn + s)*(o/Kc + 1)) - ...
        (a*mun*n*s)/(Yn*(Kan + a)*(Kn + n)^2*(Ksn + s)*(o/Kc + 1));
    fp(1,4) = (muo*o*s)/(Yo*(Kao + a)*(Ko + o)*(Kso + s)) - ...
        (a*muo*o*s)/(Yo*(Kao + a)^2*(Ko + o)*(Kso + s)) + ...
        (mun*n*s)/(Yn*(Kan + a)*(Kn + n)*(Ksn + s)*(o/Kc + 1)) - ...
        (a*mun*n*s)/(Yn*(Kan + a)^2*(Kn + n)*(Ksn + s)*(o/Kc + 1));

    fp(2,1) = (a*gamma*muo*o)/((Kao + a)*(Ko + o)*(Kso + s)) - ...
        (a*gamma*muo*o*s)/((Kao + a)*(Ko + o)*(Kso + s)^2);
    fp(2,2) = (alphao*ko)/(Kop + o) - (alphao*ko*o)/(Kop + o)^2 + ...
        (Dob*beta)/(delta*mc) + (a*gamma*muo*s)/((Kao + a)*(Ko + o)*(Kso + s)) - ...
        (a*gamma*muo*o*s)/((Kao + a)*(Ko + o)^2*(Kso + s));
    fp(2,3) = 0;
    fp(2,4) = (gamma*muo*o*s)/((Kao + a)*(Ko + o)*(Kso + s)) - ...
        (a*gamma*muo*o*s)/((Kao + a)^2*(Ko + o)*(Kso + s));

    fp(3,1) = ((a*eta*mun*n)/((Kan + a)*(Kn + n)*(Ksn + s)) - ...
        (a*eta*mun*n*s)/((Kan + a)*(Kn + n)*(Ksn + s)^2))/(o/Kc + 1);
    fp(3,2) = -((alphan*kn*n)/(Knp + n) + (a*eta*mun*n*s)/...
        ((Kan + a)*(Kn + n)*(Ksn + s)))/(Kc*(o/Kc + 1)^2);
    fp(3,3) = ((alphan*kn)/(Knp + n) - (alphan*kn*n)/(Knp + n)^2 + ...
        (a*eta*mun*s)/((Kan + a)*(Kn + n)*(Ksn + s)) - ...
        (a*eta*mun*n*s)/((Kan + a)*(Kn + n)^2*(Ksn + s)))/(o/Kc + 1) + ...
        (Dnb*beta)/(delta*mc);
    fp(3,4) = ((eta*mun*n*s)/((Kan + a)*(Kn + n)*(Ksn + s)) - ...
        (a*eta*mun*n*s)/((Kan + a)^2*(Kn + n)*(Ksn + s)))/(o/Kc + 1);

    fp(4,1) = (a*muo*o*psi)/((Kao + a)*(Ko + o)*(Kso + s)) + ...
        (a*epsilon*mun*n)/((Kan + a)*(Kn + n)*(Ksn + s)*(o/Kc + 1)) - ...
        (a*muo*o*psi*s)/((Kao + a)*(Ko + o)*(Kso + s)^2) - ...
        (a*epsilon*mun*n*s)/((Kan + a)*(Kn + n)*(Ksn + s)^2*(o/Kc + 1));
    fp(4,2) = (a*muo*psi*s)/((Kao + a)*(Ko + o)*(Kso + s)) - ...
        (a*muo*o*psi*s)/((Kao + a)*(Ko + o)^2*(Kso + s)) - ...
        (a*epsilon*mun*n*s)/(Kc*(Kan + a)*(Kn + n)*(Ksn + s)*(o/Kc + 1)^2);
    fp(4,3) = (a*epsilon*mun*s)/((Kan + a)*(Kn + n)*(Ksn + s)*(o/Kc + 1)) - ...
        (a*epsilon*mun*n*s)/((Kan + a)*(Kn + n)^2*(Ksn + s)*(o/Kc + 1));
    fp(4,4) = (Dab*beta)/(delta*mc) + (muo*o*psi*s)/((Kao + a)*(Ko + o)*(Kso + s)) + ...
        (epsilon*mun*n*s)/((Kan + a)*(Kn + n)*(Ksn + s)*(o/Kc + 1)) - ...
        (a*muo*o*psi*s)/((Kao + a)^2*(Ko + o)*(Kso + s)) - ...
        (a*epsilon*mun*n*s)/((Kan + a)^2*(Kn + n)*(Ksn + s)*(o/Kc + 1));
    
    varargout{1} = fp;
end