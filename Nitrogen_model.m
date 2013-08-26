function Nitrogen_model()

% initial values, unit: mg/l
NH4  = 0.0;
NO2  = 0;
NO3  = 5.0;
N2   = 0.0;
CH2O = 1.0;
O2   = 2.0;
X1   = 3.64/3;
X2   = 3.64/3;
X3   = 3.64/3;

% total simulation time
tfinal = 8;     % unit: day

[T, Y] = ode45(@myode, [0 tfinal],[NH4 NO2 NO3 N2 CH2O O2 X1 X2 X3]);

% plot
specs  = {'NH4', 'NO2', 'NO3', 'N2', 'CH2O', 'O2', 'X1', 'X2', 'X3'};
flag = [];
for i = 1 : 9
    if any(Y(:,i) > 1e-3)
        flag = [flag, i];
    end
end
plot(T, Y(:,flag), 'LineWidth', 2)
legend(specs(flag))
end



%--------------------------------------------------------------------------
function dydt = myode(t,y)

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
R_NH4        = 1;       % retardation factor for NH4+
R_NO2        = 1;       % retardation factor for NO2-
R_NO3        = 1;       % retardation factor for NO3-
R_N2         = 1;       % retardation factor for N2
R_CH2O       = 1;       % retardation factor for CH2O
R_O2         = 1;       % retardation factor for O2

NH4  = y(1);
NO2  = y(2);
NO3  = y(3);
N2   = y(4);
CH2O = y(5);
O2   = y(6);
X1   = y(7);    % unit: mg/l, concentration of autotrophic ammonia-oxidizing biomass
X2   = y(8);    % unit: mg/l, concentration of autotrophic nitrite-oxidizing biomass
X3   = y(9);    % unit: mg/l, the heterotrophic biomass concentration

r1 = mu_max_nit1*X1*(kb1/(kb1+X1))*(NH4/(K_NH4+NH4))*(O2/(K_O2+O2));    % unit: mg/day, substrate utilization rate by ammonium oxidation process 1
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

dydt = zeros(9,1);
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