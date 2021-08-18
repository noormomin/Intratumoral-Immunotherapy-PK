function [p,y0] = Inputs(Kd,MW,Th_L,Th_C,A,receptor,cells,NR,expand,kon_R,koff_R,kendo_R,L,C)
%INTERNAL CALL FUNCTIONS
%   %   SchmidtPerm
%   %   SchmidtVoid

%INPUT PARAMETERS
%   %   Kd          Antigen/Payload Affinity - [M]
%   %   MW          Payload Molecular Weight - [kDa]
%   %   Th_L        Payload serum half-life  - [days]
%   %   Th_C        Antigen turnover at QSS, intratumoral half life - [days]
%   %   A           Payload Uptake via Antigen Turnover - (0:no; 1:yes)
%   %   receptor    Receptor Sink Presence  - (0:none; 1:yes)
%   %   cells       Intrautmoral Receptor+ Cells - [cells/mL]
%   %   NR          Intrautmoral Receptor Density on Cells - [receptors/cell]
%   %   expand      Receptor Engaged Cell Proliferation/Expansion - [1/day]
%   %   kon_R       Payload/Receptor on rate - [1/M/s]
%   %   koff_R      Payload/Receptor off rate - [1/s]
%   %   kendo_R     Payload/Receptor endocytic rate - [1/s]
%   %   L           Payload initial concentration - [M]
%   %   C           Antigen initial concentration - [M]

%OUTPUT
%   %   p           p matrix for ODE solver
%   %   y0          y0 matrix for ODE solver 

% Initial concentrations in [M]
LC = 0;
Lesc = 2.5E-8;

Cint = 0;
LCint = 0;
R = receptor*cells*NR*10^3/(6.022E23);        % [M] - receptor concentration in tumor
LR = 0;
LCR = 0;
LRint = 0;
Reng = 0;


%ASSEMBLING VARIABLES
k1 = 1E5;                                       % [1/M/s] - on rate antigen-payload (kon,CL)
k2 = k1*Kd; %0.01;                              % [1/s] - off rate antigen-payload (koff,CL)
k6 = log(2)/(Th_C*24*60*60);                    % [1/s] - turnover of antigen
k12 = kon_R;                                    % [1/M/s] - on rate payload and receptor 
k13 = koff_R;                                   % [1/s] - off rate payload and receptor
k14 = kendo_R;                                  % [1/s] - endocytic rate of payload/receptor
k15 = R;                                        % [M] - initial receptor concentration
k16 = expand/(1*24*60*60);                      % [1/s] - receptor proliferation, logistic growth
k17 = 0.01/(60*60);                             % [1/s] - degradation of cytokine in serum measured
% k17 = 0.03/60;                                % [1/s] - kr, IL-2 receptor generation rate
% k18 = log(2)/(55*60);                         % [1/s] - 55 min half life of IL2R on surface

%Escape Mediated by Diffusion 
Mol_R = 10^(-0.31+0.43*log10(MW));   %[nm] Smilgies et al, 2015 et al (Table 2 PDF)
P = SchmidtPerm(Mol_R)*10^-2;        %[m/s] ~3.9E-9 estimate
e = SchmidtVoid(Mol_R);              %[   ] ~0.24 estimate
Rcap = 8E-6;                         %[m] capillary radius
Rk = 75E-6;                          %[m] kroghs cylinder radius
k3 = 2*P*Rcap/(Rk^2)*(1/e);         %out the tumor
k5 = e;


if isempty(Th_L)
k4 = (exp(-3.3+4.9/(1+exp((log(Mol_R)-1.4)/0.25))))/(60*60);     %[1/s] kclearance
else
k4 = log(2)/(Th_L*24*60*60);                                    % [1/s] kclearance
end

p = [k1, k2, k3, k4, k5, k6, k12, k13, k14, k15, A, k16,k17];
y0 = [L,C,LC,Lesc,Cint, LCint, LR, LCR, LRint,Reng, R];
end