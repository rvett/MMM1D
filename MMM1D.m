% ONE-DIMENSIONAL STEADY-STATE NON-ISOTHERMAL MACRO-HOMOGENEOUS TWO-PHASE
% MASTER MEA MODEL FOR PROTON EXCHANGE MEMBRANE FUEL CELLS
%
% Authors:
% Dr. Roman Vetter, roman.vetter@zhaw.ch, vetterro@ethz.ch
% Prof. Dr. Juergen O. Schumacher, juergen.schumacher@zhaw.ch
%
% Original publication:
% Free open reference implementation of a two-phase PEM fuel cell model
% Computer Physics Communications 234, 222-234 (2019)
% https://doi.org/10.1016/j.cpc.2018.07.023
%
% Copyright (c) 2017-2026 ZHAW-ICP
% Institute of Computational Physics (ICP)
% Zurich University of Applied Sciences (ZHAW)
% 8401 Winterthur
% Switzerland
%
% This software is subject to the accompanying license agreement.

function [IU, SOL] = MMM1D

% CONSTANTS
T_0 = 273.15; % [K] zero degrees Celsius
R = 8.31446; % [J/(mol*K)] universal gas constant
F = 96485.333; % [C/mol] Faraday constant
P_ref = 101325; % [Pa] reference pressure
T_ref = T_0+80; % [K] reference temperature

% ELECTROCHEMICAL PARAMETERS
BV = @(i_0,a,T,beta,eta) i_0.*a.*(exp(beta*2*F./(R*T).*eta)-exp(-(1-beta)*2*F./(R*T).*eta)); % [A/m^3] Butler-Volmer equation
A = @(E,T) exp(E/R.*(1/T_ref-1./T)); % [-] Arrhenius correction
i_0_HOR = @(T) 0.27e4*A(16e3,T); % [A/m^2] exchange current density of hydrogen oxidation reaction (HOR)
i_0_ORR = @(p_O2,T) 2.47e-4*(p_O2/P_ref).^0.54.*A(67e3,T); % [A/m^2] exchange current density of oxygen reduction reaction (ORR)
beta_HOR = 0.5; % [-] HOR symmetry factor
beta_ORR = 0.5; % [-] ORR symmetry factor
DeltaS_HOR = 0.104; % [J/(mol*K)] reaction entropy of HOR
DeltaS_ORR = -163.3; % [J/(mol*K)] reaction entropy of ORR

% WATER PROPERTIES
M_w = 18e-3; % [kg/mol] molar mass of water
V_w = M_w/0.978e3; % [m^3/mol] molar volume of liquid water (molar mass divided by mass density)
H_ec = 42e3; % [J/mol] molar enthalpy of evaporation/condensation
DeltaH = -285.83e3; % [J/mol] molar enthalpy of formation of liquid water
P_sat = @(T) exp(23.1963-3816.44./(T-46.13)); % [Pa] saturation pressure of water vapor
mu = @(T) 1e-3*exp(-3.63148+542.05./(T-144.15)); % [Pa*s] dynamic viscosity of liquid water

% GAS PROPERTIES
D = @(D_ref,eps_p,tau,s,P,T) eps_p./tau.^2.*(1-s).^3.*D_ref.*(T/T_ref).^1.5.*(P_ref./P); % [m^2/s] gas diffusivity
D_H2_ref = 1.24e-4; % [m^2/s] hydrogen diffusivity at T_ref and P_ref
D_O2_ref = 0.28e-4; % [m^2/s] oxygen diffusivity at T_ref and P_ref
D_H2O_ref_A = 1.24e-4; % [m^2/s] water vapor diffusivity in hydrogen at T_ref and P_ref
D_H2O_ref_C = 0.36e-4; % [m^2/s] water vapor diffusivity in air at T_ref and P_ref

% OPERATING CONDITIONS
U = 1.15:-0.05:0; % [V] list of cell voltages
P_A = 1.5e5; % [Pa] total pressure in anode gas channel (AGC)
P_C = 1.5e5; % [Pa] total pressure in cathode gas channel (CGC)
RH_A = 0.90; % [-] relative humidity in AGC
RH_C = 0.90; % [-] relative humidity in CGC
T_A = T_0+70; % [K] temperature of anode bipolar plate (ABP) and AGC
T_C = T_0+70; % [K] temperature of cathode bipolar plate (CBP) and CGC
alpha_H2 = 1; % [-] mole fraction of hydrogen in dry fuel gas
alpha_O2 = 0.21; % [-] mole fraction of oxygen in dry oxidant gas
s_C = 0.12; % [-] liquid water saturation at CGDL/CGC interface

x_H2O_A = RH_A*P_sat(T_A)/P_A; % [-] mole fraction of water vapor in AGC
x_H2O_C = RH_C*P_sat(T_C)/P_C; % [-] mole fraction of water vapor in CGC
x_H2_A = alpha_H2*(1-x_H2O_A); % [-] mole fraction of hydrogen in AGC
x_O2_C = alpha_O2*(1-x_H2O_C); % [-] mole fraction of oxygen in CGC

% MATERIAL PARAMETERS
L_GDL = 160e-6; % [m] GDL thickness
L_CL = 10e-6; % [m] CL thickness
L_PEM = 25e-6; % [m] PEM thickness
a_ACL = 1e7; % [1/m] electrochemical surface area density of ACL
a_CCL = 3e7; % [1/m] electrochemical surface area density of CCL
k_GDL = 1.6; % [W/(m*K)] thermal conductivity of GDL
k_CL = 0.27; % [W/(m*K)] thermal conductivity of CL
k_PEM = 0.3; % [W/(m*K)] thermal conductivity of PEM
eps_i_CL = 0.3; % [-] ionomer volume fraction in dry CL
eps_p_GDL = 0.76; % [-] porosity of GDL
eps_p_CL = 0.4; % [-] porosity of CL
kappa_GDL = 6.15e-12; % [m^2] absolute permeability of GDL
kappa_CL = 1e-13; % [m^2] absolute permeability of CL
sigma_e_GDL = 1250; % [S/m] electrical conductivity of GDL
sigma_e_CL = 350; % [S/m] electrical conductivity of CL
tau_GDL = 1.6; % [-] pore tortuosity of GDL
tau_CL = 1.6; % [-] pore tortuosity of CL
V_m = 1.02/1.97e3; % [m^3/mol] molar volume of dry membrane (equivalent weight divided by mass density)
H_ad = H_ec; % [J/mol] molar enthalpy of absorption/desorption
s_im = s_C; % [-] immobile liquid water saturation

% CONSTITUTIVE RELATIONSHIPS
iff = @(condition,a,b) condition.*a + ~condition.*b; % vectorized ternary operator
f = @(lambda) lambda*V_w./(V_m+lambda*V_w); % [-] volume fraction of water in ionomer
k_ad = @(lambda,lambda_eq,T) iff(lambda<lambda_eq,3.53e-5,1.42e-4).*f(lambda).*A(20e3,T); % [m/s] mass transfer coefficient for vapor sorption of Nafion
s_red = @(s) (s-s_im)./(1-s_im); % [-] reduced liquid water saturation
gamma_ec = @(x_H2O,x_sat,s,T) 2e6*iff(x_H2O<x_sat,5e-4*s_red(s),6e-3*(1-s_red(s))).*sqrt(R*T/(2*pi*M_w)); % [1/s] evaporation/condensation rate
sigma_p  = @(eps_i,lambda,T) eps_i.^1.5.*116.*max(0,f(lambda)-0.06).^1.5.*A(15e3,T); % [S/m] proton conductivity of Nafion
sorption = @(RH) 0.043+17.81*RH-39.85*RH.^2+36.0*RH.^3; % [-] vapor sorption isotherm of Nafion
D_lambda = @(eps_i,lambda,T) eps_i.^1.5.*(3.842*lambda.^3-32.03*lambda.^2+67.74*lambda)./(lambda.^3-2.115*lambda.^2-33.013*lambda+103.37)*1e-10.*A(20e3,T); % [m^2/s] diffusion coefficient of water in Nafion
xi = @(lambda) 2.5*lambda/22; % [-] electro-osmotic drag coefficient of Nafion
dpds = @(s) 0.00011*44.02*exp(-44.02*(s-0.496))+278.3*8.103*exp(8.103*(s-0.496)); % [Pa] derivative of capillary pressure-saturation relationship of GDL and CL
D_s = @(kappa,s,T) kappa.*(1e-6+s_red(s).^3)./mu(T).*dpds(s); % [m^2/s] liquid water transport coefficient

% PER-LAYER PROPERTIES (columns: AGDL, ACL, PEM, CCL, CGDL)
prop = [L_GDL       L_CL        L_PEM L_CL        L_GDL;        % [m] thickness
        0           eps_i_CL    1     eps_i_CL    0;            % [-] ionomer volume fraction
        eps_p_GDL   eps_p_CL    0     eps_p_CL    eps_p_GDL;    % [-] porosity
        tau_GDL     tau_CL      1     tau_CL      tau_GDL;      % [-] pore tortuosity
        sigma_e_GDL sigma_e_CL  0     sigma_e_CL  sigma_e_GDL;  % [S/m] electrical conductivity
        k_GDL       k_CL        k_PEM k_CL        k_GDL;        % [W/(m*K)] thermal conductivity
        kappa_GDL   kappa_CL    0     kappa_CL    kappa_GDL;    % [m^2] absolute permeability
        P_A         P_A         P_ref P_C         P_C;          % [Pa] gas pressure
        D_H2O_ref_A D_H2O_ref_A 0     D_H2O_ref_C D_H2O_ref_C]; % [m^2/s] water vapor diffusivity prefactor

Nd = size(prop,2); % number of domains (layers)
idx = num2cell(1:Nd);
[AGDL, ACL, PEM, CCL, CGDL] = idx{:};

% DEPENDENT QUANTITIES (each has a potential-like variable and a flux in y; rows interleaved)
quantities = {'phi_e', 'phi_p', 'T', 'lambda', 'x_H2O', 'x_H2', 'x_O2', 's'};
Nq = numel(quantities); % number of quantities
q = cell2struct(num2cell(1:Nq    ), quantities, 2); % quantity index
v = cell2struct(num2cell(1:2:2*Nq), quantities, 2); % index of the field variable

% ACTIVE LAYERS OF EACH QUANTITY
domains = false(Nq, Nd);
domains(q.phi_e , [AGDL ACL     CCL CGDL]) = true;
domains(q.phi_p , [     ACL PEM CCL     ]) = true;
domains(q.T     , [AGDL ACL PEM CCL CGDL]) = true;
domains(q.lambda, [     ACL PEM CCL     ]) = true;
domains(q.x_H2O , [AGDL ACL     CCL CGDL]) = true;
domains(q.x_H2  , [AGDL ACL             ]) = true;
domains(q.x_O2  , [             CCL CGDL]) = true;
domains(q.s     , [             CCL CGDL]) = true;
active = repelem(domains, 2, 1); % active fields and fluxes (2*Nq x Nd)

% INITIAL GUESS FOR FIRST CELL VOLTAGE
y0 = zeros(2*Nq, Nd); % all fluxes and phi_p = 0
y0(v.phi_e , [             CCL CGDL]) = U(1);
y0(v.T     , [AGDL ACL PEM CCL CGDL]) = (T_A+T_C)/2;
y0(v.lambda, [     ACL PEM CCL     ]) = sorption(1);
y0(v.x_H2O , [AGDL ACL             ]) = x_H2O_A;
y0(v.x_H2O , [             CCL CGDL]) = x_H2O_C;
y0(v.x_H2  , [AGDL ACL             ]) = x_H2_A;
y0(v.x_O2  , [             CCL CGDL]) = x_O2_C;
y0(v.s     , [             CCL CGDL]) = s_C;

% INITIAL MESH
Lsum = [0 cumsum(prop(1,:))];
Nx = 10; % number of initial mesh intervals per domain
x = reshape((0:Nx)'/Nx*prop(1,:) + Lsum(1:Nd), 1, []); % interface nodes appear twice

% SOLVER PREPARATION
sol = bvpinit(x, @(~, d) y0(:, d));
options = bvpset('Vectorized', 'on', 'NMax', 100*Nd, 'RelTol', 1e-4, 'AbsTol', 1e-6);

% CELL VOLTAGE SWEEP
SOL = cell(size(U));
for m = 1:numel(U)
    sol = bvp4c(@odefun, @(ya,yb) bcfun(ya, yb, m), sol, options);
    SOL{m} = sol;
end
I = cellfun(@(s) s.y(v.phi_e+1,end), SOL)/1e4; % current density in [A/cm^2]
IU = [I(:) U(:)];

% POSTPROCESSING
Nref = 2; % number of refinements for smoother curve plotting
offset = 1e-10; % a small offset (relative to interval width) to evaluate solutions just inside interval ends
pos = linspace(offset, 1-offset, 2^Nref+1)'; % relative positions within each mesh interval
for m = 1:numel(U)
    x = interp1(SOL{m}.x, reshape(find(diff(SOL{m}.x) > 0) + pos, 1, [])); % refined points of all non-degenerate mesh intervals
    SOL{m}.y = deval(SOL{m}, x);
    SOL{m}.y(~active(:, discretize(x, Lsum))) = NaN; % NaN on inactive domains
    SOL{m}.x = x;
end

% PLOT ALL FIELDS & FLUXES
% plotting info per quantity: {field label, flux label, flux unit scale}
info.phi_e  = {'\phi_e [V]', 'j_e [A/cm^2]'            , 1e-4};
info.phi_p  = {'\phi_p [V]', 'j_p [A/cm^2]'            , 1e-4};
info.T      = {'T [K]'     , 'j_T [W/cm^2]'            , 1e-4};
info.lambda = {'\lambda'   , 'j_\lambda [\mumol/cm^2s]', 1e2 };
info.x_H2O  = {'x_{H_2O}'  , 'j_{H_2O} [\mumol/cm^2s]' , 1e2 };
info.x_H2   = {'x_{H_2}'   , 'j_{H_2} [\mumol/cm^2s]'  , 1e2 };
info.x_O2   = {'x_{O_2}'   , 'j_{O_2} [\mumol/cm^2s]'  , 1e2 };
info.s      = {'s'         , 'j_s [\mumol/cm^2s]'      , 1e2 };

fig_names = {'Fields', 'Fluxes'};
Nr = ceil(sqrt(Nq)); Nc = ceil(Nq/Nr); % number of rows and columns in subplot grid
c = jet(numel(U)); % colors
for p = 1:numel(fig_names)
    figure('Name', fig_names{p})
    for n = 1:Nq
        name = quantities{n};
        scale = [1 info.(name){3}];
        subplot(Nr, Nc, n)
        box on
        hold on
        for m = 1:numel(U)
            plot(SOL{m}.x*1e6, SOL{m}.y(v.(name)+p-1,:)*scale(p), 'Color', c(m,:), 'DisplayName', [num2str(U(m)) ' V'])
        end
        d = find(domains(n,:)); % domains in which this quantity is active
        xlim(Lsum([d(1) d(end)+1])*1e6)
        xlabel('x [\mum]')
        ylabel(info.(name){p})
        xline(Lsum(2:end-1)*1e6, 'k', 'HandleVisibility', 'off');
    end
end

% PLOT POLARIZATION CURVE
figure('Name', 'Polarization curve')
box on
xlabel('Current density [A/cm^2]')

yyaxis left
plot(I, U)
ylabel('Cell voltage [V]')

yyaxis right
plot(I, I.*U)
ylabel('Power density [W/cm^2]')

% SYSTEM OF 1ST-ORDER DIFFERENTIAL EQUATIONS
function dydx = odefun(x, y, d)
    
    pr = num2cell(prop(:,d));
    [L, eps_i, eps_p, tau, sigma_e, k, kappa, P, D_H2O_ref] = pr{:}; % properties of this domain
    
    phi_e  = y(v.phi_e  ,:);  j_e = y(v.phi_e+1,:);
    phi_p  = y(v.phi_p  ,:);  j_p = y(v.phi_p+1,:);
    T      = y(v.T      ,:);
    lambda = y(v.lambda ,:);
    x_H2O  = y(v.x_H2O  ,:);
    x_H2   = y(v.x_H2   ,:);
    x_O2   = y(v.x_O2   ,:);
    s      = y(v.s      ,:);
    
    C = P./(R*T); % total interstitial gas concentration
    x_sat = P_sat(T)./P; % saturation water vapor mole fraction
    
    % FLUX LAWS: j_q = -K_q*grad(q) (plus electro-osmotic drag for lambda)
    K = zeros(Nq,numel(x));
    K(q.phi_e ,:) = sigma_e;
    K(q.phi_p ,:) = sigma_p(eps_i,lambda,T);
    K(q.T     ,:) = k;
    K(q.lambda,:) = D_lambda(eps_i,lambda,T)/V_m;
    CD = @(D_ref) C.*D(D_ref,eps_p,tau,s,P,T); % gas flux coefficient
    K(q.x_H2O ,:) = CD(D_H2O_ref);
    K(q.x_H2  ,:) = CD(D_H2_ref);
    K(q.x_O2  ,:) = CD(D_O2_ref);
    K(q.s     ,:) = D_s(kappa,s,T)/V_w;
    K(~domains(:,d),:) = Inf; % inactive quantities: zero gradient
    
    j = y(2:2:end,:); % fluxes
    j(q.lambda,:) = j(q.lambda,:)-xi(lambda)/F.*j_p; % electro-osmotic drag
    dydx = zeros(size(y));
    dydx(1:2:end,:) = -j./K; % field gradients
    
    % ELECTROCHEMICAL REACTION
    eta = zeros(size(x)); % overpotential
    i = zeros(size(x)); % volumetric current density
    nu = zeros(Nq,1); % stoichiometry per unit volumetric current density
    DeltaS = 0;
    switch d
        case ACL % hydrogen oxidation reaction
            DeltaS = DeltaS_HOR;
            eta = phi_e-phi_p + T*DeltaS/(2*F) + R*T/(2*F).*log(x_H2.*P/P_ref);
            i = BV(i_0_HOR(T),a_ACL,T,beta_HOR,eta);
            nu([q.phi_e q.phi_p q.x_H2]) = [-1 1 -1/(2*F)]';
        case CCL % oxygen reduction reaction
            DeltaS = DeltaS_ORR;
            eta = -(phi_e-phi_p) - (DeltaH-T*DeltaS)/(2*F) + R*T/(4*F).*log(x_O2.*P/P_ref);
            i = BV(i_0_ORR(x_O2.*P,T),a_CCL,T,beta_ORR,eta);
            nu([q.phi_e q.phi_p q.lambda q.x_O2]) = [1 -1 1/(2*F) -1/(4*F)]';
    end
    
    % CONSERVATION LAWS: div(j_q) = S_q
    S = nu*i; % reaction sources
    S_Te = -j_e.*dydx(v.phi_e,:); % Joule's 1st law
    S_Tp = -j_p.*dydx(v.phi_p,:);
    S_Tr = i.*(eta-T*DeltaS/(2*F)); % reaction heat
    lambda_eq = sorption(x_H2O./x_sat); % equilibrium water content of ionomer
    S_ad = (domains(q.lambda,d)&&domains(q.x_H2O,d))*k_ad(lambda,lambda_eq,T)/(L*V_m).*(lambda_eq-lambda); % absorption/desorption
    S_ec = (domains(q.s,d)&&domains(q.x_H2O,d))*gamma_ec(x_H2O,x_sat,s,T).*C.*(x_H2O-x_sat); % evaporation/condensation
    
    S(q.T       ,:) = S_Te + S_Tp + S_Tr + H_ad.*S_ad + H_ec.*S_ec; % Joule heat + reaction heat + phase-change heat
    S(q.lambda  ,:) = S(q.lambda,:) + S_ad;
    S(q.x_H2O   ,:) = -S_ec - S_ad;
    S(q.s       ,:) = S_ec;
    dydx(2:2:end,:) = S; % sources
end

% BOUNDARY & INTERFACE CONDITIONS
function res = bcfun(ya, yb, m)
    
    jump = ya(:,2:end) - yb(:,1:end-1);      % right minus left value at each interface
    isflux = repmat([false; true], Nq, 1);   % rows that are fluxes
    al = active(:,1:end-1);                  % active left of each interface?
    ar = active(:,2:end);                    % active right of each interface?
    link = (al & ar) | (isflux & (al | ar)); % variables tied across an interface
    
    res = [ya(~active);                % inactive variables pinned to zero
           jump(link);                 % continuity (= zero flux where a quantity ends or starts)
           ya(v.phi_e, 1  );           % phi_e = 0 (electrical ground)
           ya(v.T    , 1  ) - T_A;     % ABP temperature
           ya(v.x_H2O, 1  ) - x_H2O_A; % AGC vapor mole fraction
           ya(v.x_H2 , 1  ) - x_H2_A;  % AGC hydrogen mole fraction
           yb(v.phi_e, end) - U(m);    % cell voltage
           yb(v.T    , end) - T_C;     % CBP temperature
           yb(v.x_H2O, end) - x_H2O_C; % CGC vapor mole fraction
           yb(v.x_O2 , end) - x_O2_C;  % CGC oxygen mole fraction
           yb(v.s    , end) - s_C];    % CGC liquid water saturation
end

end
