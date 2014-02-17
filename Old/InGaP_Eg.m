function [Eg, ni, muw, mue, mui, mub]=InGaP_Eg(GaInP_comp, Ne_a, Nb_d, T)



GaP_Eg0_Gam = 2.886;
GaP_alp_Gam = 0.1081;
GaP_Beta_Gam = 164;
GaP_Eg0_X = 2.35;
GaP_alp_X = 0.5771e-3;
GaP_Beta_X = 372;
GaP_Eg0_L = 2.72;
GaP_alp_L = 0.5771e-3;
GaP_Beta_L = 372;

InP_Eg0_Gam = 1.4236;
InP_alp_Gam = 0.363e-3;
InP_Beta_Gam = 162;
InP_Eg0_X = 2.384;
InP_alp_X = 0.37e-3;
InP_Eg0_L = 2.014;
InP_alp_L = 0.363e-3;
InP_Beta_L = 162;

GaInP_C_Gam = 0.65;
GaInP_C_X = 0.20;
GaInP_C_L = 1.03;


GaP_Eg_T_Gam = GaP_Eg0_Gam + GaP_alp_Gam * (1 - coth(GaP_Beta_Gam/T));
GaP_Eg_T_X = GaP_Eg0_X - GaP_alp_X * T^2 / (T + GaP_Beta_X);
GaP_Eg_T_L = GaP_Eg0_L - GaP_alp_L * T^2 / (T + GaP_Beta_L);

InP_Eg_T_Gam = InP_Eg0_Gam - InP_alp_Gam * T^2 / (T + InP_Beta_Gam);
InP_Eg_T_X = InP_Eg0_X - InP_alp_X * T;
InP_Eg_T_L = InP_Eg0_L - InP_alp_L * T^2 / (T + InP_Beta_L);



    GaInP_Eg_Gam = (1-GaInP_comp) * GaP_Eg_T_Gam + GaInP_comp* InP_Eg_T_Gam ...
        - GaInP_comp * (1-GaInP_comp) * GaInP_C_Gam
    GaInP_Eg_X = (1-GaInP_comp) .* GaP_Eg_T_X + GaInP_comp .* InP_Eg_T_X ...
        - GaInP_comp * (1-GaInP_comp) * GaInP_C_X
    GaInP_Eg_L = (1-GaInP_comp) .* GaP_Eg_T_L + GaInP_comp .* InP_Eg_T_L ...
        - GaInP_comp * (1-GaInP_comp) * GaInP_C_L
    GaInP_Eg_compare=[GaInP_Eg_Gam, GaInP_Eg_X, GaInP_Eg_L]
    GaInP_Eg_T = min(GaInP_Eg_compare);
   

Eg=GaInP_Eg_T

kb=1.381e-23;
q=1.602000000000000e-19;
Vt=kb*T/q;
erel_e=11.1;%erel_GaAs;
erel_i=11.1;%erel_GaAs; 
erel_b=11.1;%erel_GaAs;

Nc=6.235E17*(T/300)^(3/2);%LM to GaAs
Nv=1.457E19*(T/300)^(3/2);%LM to GaAs
ni=((Nc*Nv)^.5)*exp((-Eg)./(2*Vt));


Ncr=1.6e24*(.066/(1.4*erel_b))^3;
Maw=2E18; %maj carrier conc in window
Ma=Ne_a; %maj carrier conc in emitter
Mi=ni;
Md=Nb_d; %maj carrier conc in base
Msub=1E18;

c=1/3;
mumaxn=-6933*c^3+3581*c^2+76*c+9157;
muminn=-2201*c^3+5362*c^2-5092*c+1942;
Nrefn=-4.36e16*c^3+1.43e17*c^2-1.80e17*c+8.49e16;
chi=.506*c^3-.262*c^2+.132*c+.365;
Tn=T/300;

mumaxp=150; %[cm2/Vs]
Nrefp=6.25e17; %[cm-3]

muw=muminn*Tn^(-.57)+((mumaxn-muminn)*Tn^(-2.3))/(1+(Maw/(Nrefn*Tn^(-2.4)))^(chi*Tn^(-.146))); %InGaP window  (minority electrons in ptype)
mue=muminn*Tn^(-.57)+((mumaxn-muminn)*Tn^(-2.3))/(1+(Maw/(Nrefn*Tn^(-2.4)))^(chi*Tn^(-.146))); %InGaP window  (minority electrons in ptype)
mui=muminn*Tn.^(-.57)+((mumaxn-muminn)*Tn^(-2.3))./(1+(Mi./(Nrefn*Tn^(-2.4))).^(chi*Tn^(-.146))); %pGaAs i-region
mub=mumaxp*(Tn.^(-2.3))./(1+Md./(Nrefp*Tn^3.8)) %nGaAs base (minority holes in n type)
musub=1.5E4;%mumaxp*(Tn^(-2.3))/(1+Msub/(Nrefp*Tn^3.8));
