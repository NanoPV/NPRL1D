function [Eg, ni, mu, meeff,mheff erel]=nInGaP_Eg(GaInP_comp, N, T)
%In xIn
%% Constants
kb=1.381e-23;%J/K
k=8.6173324e-5; %[eV/K] (Wikipedia)
q=1.602176565e-19;%C
Vt=kb*T/q; %V
h=6.62606957e-34; %Js
hbar=h/(2*pi); %[J.s] (obvs)
m0=9.11E-31; %kg
epsilon_0=8.854187817620e-12; %[F/m]
Zd=1; %Absolute value of ion charge of acceptor
Za=1; %Absolute value of ion charge of donor
m_to_cm=1e-6; % Conversion from m^-3 to cm^-3
cm_to_m=1e6; % Conversion fromc cm^-3 to m^-3
cm2_to_m2=1e-4; % Conversion from cm^2 to m^2
m2_to_cm2=1e4; % Conversion from m^2 to cm^2
dyn_to_N=1e-5; % Conversion from dyne to Newton
gd=2; %Dengeneracy of electrons in conduction band
ga=4; %Degeneracy of heavy holes and light holes in valence band
Tn=T/300; %Normalized T
%% Bowing Parameters for effective mass
F0=-2.04;
F1=-1.31;
Fbow=0.78;

Ep0=31.4;
Ep1=20.7;
Epbow=0;

dSO0=0.08;
dSO1=0.108;
dSObow=0;

mh0=0.6;
mh1=0.64;;
%% Bandgap Calculation
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
Ec=Eg/2; %[eV]
Ev=-Ec; %[eV]
Ed=Ec;
Ea=Ec;


erel=11.1+1.4*GaInP_comp;%relative dielectric permittivity
%% Effective Mass Calculation
Fcomp=F0+(F1-F0)*GaInP_comp+Fbow*GaInP_comp*(1-GaInP_comp);
Epcomp=Ep0+(Ep1-Ep0)*GaInP_comp+Epbow*GaInP_comp*(1-GaInP_comp);
dSOcomp=dSO0+(dSO1-dSO0)*GaInP_comp+Epbow*GaInP_comp*(1-GaInP_comp);
meeff=1/((1+2*Fcomp)+Epcomp*(Eg+2*dSOcomp/3)/(Eg*(Eg+dSOcomp)));
mheff=mh0+(mh1-mh0)*GaInP_comp;

%% Density of States, Ni, and Fermi Level
Nc=2*(2*pi*meeff*m0*kb*T/(h^2))^(3/2)*1E-6;
Nv=2*(2*pi*mheff*m0*kb*T/(h^2))^(3/2)*1E-6;
ni=((Nc*Nv)^.5)*exp((-Eg)./(2*Vt));

Efn=fsolve(@(Efn) fermi_charge_neutrality(Efn,Ec,k,T,Ev,Nc,Nv,0,ga,Ea,N,gd,Ed),Ec-1e-8,optimset('Display','off')) %[eV]

%% Mobility
mumaxp=150; %[cm2/Vs]
Nrefp=6.25e17; %[cm-3]

mu=mumaxp*(Tn.^(-2.3))./(1+N/(Nrefp*Tn^3.8)) %nGaAs base (minority holes in n type)

