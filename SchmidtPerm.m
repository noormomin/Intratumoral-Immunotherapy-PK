%   Written by Noor Momin
%   Last Updated: November 23, 2019 by Noor Momin
%%   References
%   (1) Schmidt and Wittrup, Mol Cancer Ther. 2009
%   (2) Nugent and Jain et al 1983 
%   (3) Paine and Scherr et al 1975
%%  Input and Outputs
%   Input: Mol_R - molecular radius in [nm]
%   Output: P - tumor capillary permeability  [cm/s]
%   Assumptions: Fit parameters (ref 1) for two-pore model of the capillary wall
%%
function P = SchmidtPerm(Mol_R)
%Paine and Scherr look up values
PS = [0.60	11.14580;
      0.62	12.87453;
      0.64	14.98751;
      0.66	17.59862;
      0.68	20.86550;
      0.70	25.01092;
      0.72	30.35733;
      0.74	37.38467;
      0.76	46.83132;
      0.78	59.87967;
      0.80	78.51925;
      0.82	106.31670;
      0.84	150.22684;
      0.86	225.51931;
      0.88	372.41035;
      0.90	737.25652]; 

%Two-Pore Model of Capillary Wall Fit Parameters (ref 1)
Rsmall = 4.5;       % small capillary radius [nm]
fasmall = 17.6;     % small capillary fractional area to thickness ratio [cm^-1]
Rlarge = 500;       % large capillary radius [nm]
falarge = 0.65;     % large capillary fractional area to thickness ratio [cm^-1]
Dfree = 3E-6/Mol_R; % diffusivity of molecule in solution estimate [cm^2/s]

%For the Small Capillary
L = Mol_R/Rsmall;   %ratio of molecular radius to pore radius

    %Conditional Dpore/Dfree (DpDf) and partition coefficient (sig)
    if L<0.6 
        DpDfsmall = (1 - 2.105*L + 2.0865*L^3 - 1.7068*L^5 + 0.72603*L^6)/(1-0.75857*L^5);
        sigsmall = (1-L)^2;
    elseif 0.6<=L && L<= 0.9
        [minD, minI] = min(abs(L-PS(:,1)));
        DpDfsmall = (1-L)^2/PS(minI,2);
        sigsmall = (1-L)^2;
    elseif 0.9<L
        DpDfsmall = 0;
        sigsmall = 0;
    end
    
%For the Large Capillary
L = Mol_R/Rlarge;   %ratio of molecular radius to pore radius

    %Conditional Dpore/Dfree (DpDf) and partition coefficient (sig)
    if L < 0.6 
        DpDflarge = (1 - 2.105*L + 2.0865*L^3 - 1.7068*L^5 + 0.72603*L^6)/(1-0.75857*L^5);
        siglarge = (1-L)^2;
    elseif 0.6<=L && L<= 0.9
        [minD, minI] = min(abs(L-PS(:,1)));
        DpDflarge = (1-L)^2/PS(minI,2);
        siglarge = (1-L)^2;
    elseif 0.9<L
        DpDflarge = 0;
        siglarge = 0;
    end

%Pore Permeability [cm^2/s]
Plarge = Dfree*DpDflarge*siglarge*falarge;
Psmall = Dfree*DpDfsmall*sigsmall*fasmall;

%Total Permeability 
P = Psmall + Plarge;

end