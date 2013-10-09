function Nitrogen_model3()
% parameters and model equations refer to
% Porporato et al. Hydrologic controls on soil carbon and nitrogen cycles.
% I. Modeling scheme; II. A case study
% AWR, 2003

time = (0:5000)';
dt = 1;
y0 = [0.11, 1500, 8500, 50, 37.5, 0.002, 0.3];
y = zeros(length(time), length(y0));
DECl   = zeros(length(time), 1);
MIN    = zeros(length(time), 1);
UP_NO3 = zeros(length(time), 1);
LE_NO3 = zeros(length(time), 1);
LE_NH4 = zeros(length(time), 1);
UP_NH4 = zeros(length(time), 1);
y(1, :) = y0;
for n = 1 : length(time) - 1
    [dCNdt, DECl(n), MIN(n), UP_NO3(n), LE_NO3(n), LE_NH4(n), UP_NH4(n)] = myode(time(n), y(n,:));
    y(n+1, :) = y(n, :) + dt * dCNdt;
end
legends = {'soil moisture','C_l (gC m^{-3})','C_h (gC m^{-3})','C_b (gC m^{-3})','N_l (gN m^{-3})','NH_4^+ (gN m^{-3})','NO_3^- (gN m^{-3})'};
for i = 1 : length(y0)
    figure;
    plot(time, y(:,i));
    title(legends{i});
    xlabel('t (d)');
end
figure; plot(time, DECl); title('Litter decomposition (gC m^{-3}d^{-1})'); xlabel('t (d)');
figure; plot(time, MIN); title('Net Mineralization (gN m^{-3}d^{-1})'); xlabel('t (d)')
figure; plot(time, UP_NO3); title('NO_3^- uptake (gN m^{-3}d^{-1})'); xlabel('t (d)')
figure; plot(time, LE_NO3); title('NO_3^- leaching (gN m^{-3}d^{-1})'); xlabel('t (d)')
figure; plot(time, LE_NH4); title('LE NH4');
figure; plot(time, UP_NH4); title('UP NH4');


end

function [dCNdt, DECl, MIN, UP_NO3, LE_NO3, LE_NH4, UP_NH4] = myode(t, CN)

date = datestr(t + datenum('01-Jan-1998'));
month = date(4:6);
switch month
    case 'Jan'
        ADD = 0.002;
    case 'Feb'
        ADD = 0.002;
    case 'Mar'
        ADD = 0.00227;
    case 'Apr'
        ADD = 0.00227;
    case 'May'
        ADD = 0.00255;
    case 'Jun'
        ADD = 0.00271;
    case 'Jul'
        ADD = 0.00282;
    case 'Aug'
        ADD = 0.006376;
    case 'Sep'
        ADD = 0.00473;
    case 'Oct'
        ADD = 0.00282;
    case 'Nov'
        ADD = 0.002;
    case 'Dec'
        ADD = 0.002;
    otherwise
end
ADD = ADD * 100 * 5.2561;     % convert from Mg C ha^-1 day^-1 to g C m^-2 day^-1
                              % data from Christiane AWR 2012 paper figure
                              % 4 and make the mean to table 1 on P.
                              % D'Odorico AWR 2003 paper

s  = CN(1);     % soil moisture
Cl = CN(2);     % carbon concentration in the litter pool
Ch = CN(3);     % carbon concentration in the humus pool
Cb = CN(4);     % carbon concentration in the biomass pool
Nl = CN(5);     % organic nitrogen concentration in the litter pool
% Nh = CN(6);     % organic nitrogen concentration in the humus pool
% Nb = CN(7);     % organic nitrogen concentration in the biomass pool
NH4 = CN(6);    % ammonium concentration in the soil
NO3 = CN(7);    % nitrate concentration in the soil

C_over_N_add = 58;
C_over_N_b = 11.5;
C_over_N_h = 22;
C_over_N_l = Cl / Nl;

% soil moisture part

sstar = 0.17;           % dimensionless, soil moisture when plants close stomata
sW    = 0.065;          % dimensionless, soil moisture at wilting point
sh    = 0.02;           % dimensionless, hygroscopic soil moisture
EW    = 0.01/100;       % unit: cm day^-1 --> m day^-1, evaporation rate at wilting point
EMAX  = 4.5/1000;       % unit: mm day^-1 --> m day^-1, maximum transpiration rate
if sh < s && s <= sW
    Es = EW * (s - sh) / (sW - sh);
else
    if s <= sstar
        Es = EW + (EMAX - EW) * (s - sW) / (sstar - sW);
    else
        if s <= 1
            Es = EMAX;
        end
    end
end

TMAX = 0.13/100;         % unit: cm day^-1 --> m day^-1, maximum evaporation rate
if 0 <= s && s <= sW
    Ts = 0;
else
    if s <= sstar
        Ts = (s - sW) / (sstar - sW) * TMAX;
    else
        if s <= 1
            Ts = TMAX;
        end
    end
end

n   = 0.4;          % Dimensionless, porosity
Zr  = 0.8;          % unit: m, soil depth
sfc = 0.3;          % dimensionless
Ks  = 1.1;          % unit: m day^-1, saturated hydraulic conductivity
% b is the pore size distribution index
%                 b       c     beta
% Sand          4.05    11.1    12.1
% Loamy sand    4.38    11.7    12.7
% Sandy Load    4.90    12.8    13.8
% Loam          5.39    13.8    14.8
% Clay          11.4    25.8    26.8
b = 4;        % used by the paper author
beta = 2 * b + 4;
% vertical percolation with unit gradient
% Ks is the saturated hydraulic conductivity
if s > sfc
    Ls = Ks / (exp(beta*(1-sfc))-1) * (exp(beta*(s-sfc))-1);    % sfc < s <=1
else
    Ls = 0;
end
% soil moisture dynamics, Ist is the rate of infiltration from rainfall
Ist = rainfall();
Ist = min([Ist, n*Zr*(1-s)]);
dsdt = (Ist - Es - Ts - Ls) / (n * Zr);

% C part

kd = 8.5e-3;    % unit: d^-1, rate of carbon return to litter pool due to death of microbial biomass
BD = kd * Cb;

if s <= sfc         % sfc is soil moisture at field capacity
    fd = s / sfc;   % nondimensional factor describes the effects of soil
                    % moisture on decomposition
else
    fd = sfc / s;
end

kiNH4 = 1;      % unit: m^3 d^-1 g C^-1, immobilization rate for NH4
kiNO3 = 1;      % unit: m^3 d^-1 g C^-1, immobilization rate for NO3
kh    = 2.5e-6; % unit: m^3 d^-1 g N^-1
rr = 0.6;       % dimensionless, fraction of decomposing C lost to respiration
rh = min(0.25, C_over_N_h / C_over_N_l);
kl = 6.5e-5;    % unit: m^3 d^-1 g C^-1

eq18_curly_brackets = kh*Ch*(1/C_over_N_h-(1-rr)/C_over_N_b)+kl*Cl*...
    (1/C_over_N_l-rh/C_over_N_h-(1-rh-rr)/C_over_N_b);

if eq18_curly_brackets > 0
    phi = 1;
else
    IMM_max = (kiNH4 * NH4 + kiNO3 * NO3) * fd * Cb;
    if - fd * Cb * eq18_curly_brackets < IMM_max
        phi = 1;
    else
        phi = - (kiNH4 * NH4 + kiNO3 * NO3) / eq18_curly_brackets;
    end
end

PHI = phi * fd * Cb * eq18_curly_brackets;

if PHI > 0
    MIN = PHI;
    IMM = 0;
    IMM_NH4 = IMM;
    IMM_NO3 = IMM;
else
    MIN = 0;
    IMM = - PHI;
    IMM_NH4 = (kiNH4 * NH4) / (kiNH4 * NH4 + kiNO3 * NO3) * IMM;
    IMM_NO3 = (kiNO3 * NO3) / (kiNH4 * NH4 + kiNO3 * NO3) * IMM;
end

DECl = phi * fd * kl * Cb * Cl;     % carbon output due to microbial dceomposition

dCldt = ADD + BD - DECl;
dNldt = ADD/C_over_N_add + BD/C_over_N_b - DECl/C_over_N_l;

DECh = phi * fd * kh * Cb * Ch;

dChdt = rh * DECl - DECh;
% dNhdt = dChdt / C_over_N_h;

dCbdt = (1 - rh - rr) * DECl + (1 - rr) * DECh - BD;
% dNbdt = (1 - rh * C_over_N_l / C_over_N_h) * DECl / C_over_N_l + ...
%     DECh / C_over_N_h - BD / C_over_N_b - PHI;

a_NH4 = 0.05;       % dimensionless
a_NO3 = 1;          % dimensionless
LE_NH4 = a_NH4 * Ls * NH4 / (s * n * Zr);
LE_NO3 = a_NO3 * Ls * NO3 / (s * n * Zr);

UPp_NH4 = a_NH4 * Ts * NH4 / (s * n * Zr);
UPp_NO3 = a_NO3 * Ts * NO3 / (s * n * Zr);

DEM_NH4 = 0.2;      % unit: g N m^-3 d^-1
DEM_NO3 = 0.5;      % unit: g N m^-3 d^-1
F = 0.1;            % unit: m d^-1, rescaled diffusion coefficient
d = 3;              % dimensionless, nonlinear dependence of the diffusion process on soil moisture
ku = a_NH4 / (s * n * Zr) * F * s^d;
if DEM_NH4 - UPp_NH4 < 0
    UPa_NH4 = 0;
else
    if DEM_NH4 - UPp_NH4 < ku * NH4
        UPa_NH4 = DEM_NH4 - UPp_NH4;
    else
        UPa_NH4 = ku * NH4;
    end
end
ku = a_NO3 / (s * n * Zr) * F * s^d;
if DEM_NO3 - UPp_NO3 < 0
    UPa_NO3 = 0;
else
    if DEM_NO3 - UPp_NO3 < ku * NO3
        UPa_NO3 = DEM_NO3 - UPp_NO3;
    else
        UPa_NO3 = ku * NO3;
    end
end

UP_NH4 = UPp_NH4 + UPa_NH4;
UP_NO3 = UPp_NO3 + UPa_NO3;

% fn(s) accounts for the soil moisture effects on nitrification
if s <= sfc
    fn = s / sfc;
else
    fn = (1 - s) / (1 - sfc);
end
kn = 0.006;       % unit: m^3 d^-1 g C^-1, the rate of nitrification
NIT = fn * kn * Cb * NH4;

dNH4dt = MIN - IMM_NH4 - NIT - LE_NH4 - UP_NH4;
dNO3dt = NIT - IMM_NO3 - LE_NO3 - UP_NO3;
% dNH4dt = 0;
% dNO3dt = 0;
% kn_NH4 = (MIN - IMM_NH4 - LE_NH4 - UP_NH4) / (fn * Cb * NH4);
% kn_NO3 = (IMM_NO3 + LE_NO3 + UP_NO3) / (fn * Cb * NH4);

%
dCNdt    = zeros(size(CN));
dCNdt(1) = dsdt;
dCNdt(2) = dCldt;
dCNdt(3) = dChdt;
dCNdt(4) = dCbdt;
dCNdt(5) = dNldt;
% dCNdt(6) = dNhdt;
% dCNdt(7) = dNbdt;
dCNdt(6) = dNH4dt;
dCNdt(7) = dNO3dt;

end