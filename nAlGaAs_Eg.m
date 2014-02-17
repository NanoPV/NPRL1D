function [Eg, ni, mu, erel]=nAlGaAs_Eg(AlGaAs_comp, N, T)
%In xAs
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
%% Bowing Parameters for effective mass
F0=-2.04;
F1=-1.94;
Fbow=0;

Ep0=31.4;
Ep1=28.8;
Epbow=0;

dSO0=0.08;
dSO1=0.341;
dSObow=0;

mh0=0.6;
mh1=0.082;
%% Bandgap Calculation
AlAs_Eg0_Gam = 3.099; %[eV]
AlAs_alp_Gam = 0.885e-3; %[eV/K]
AlAs_Beta_Gam = 530; %[K]
AlAs_Eg0_X = 2.24; %[eV]
AlAs_alp_X = 0.7e-3; %[eV/K]
AlAs_Beta_X = 530; %[K]
AlAs_Eg0_L = 2.46; %[eV]
AlAs_alp_L = 0.605e-3; %[eV/K]
AlAs_Beta_L = 204; %[K]
AlAs_Lat_Const = 5.6611 + 2.9e-5 * (Temp-300); %[Ang]

GaAs_Eg0_Gam = 1.519;
GaAs_alp_Gam = 0.5405e-3;
GaAs_Beta_Gam = 204;
GaAs_Eg0_X = 1.981;
GaAs_alp_X = 0.46e-3;
GaAs_Beta_X = 204;
GaAs_Eg0_L = 1.815;
GaAs_alp_L = 0.605e-3;
GaAs_Beta_L = 204;

GaAsP_C_Gam = 0.19;
GaAsP_C_X = 0.24;
GaAsP_C_L = 0.16;


GaP_Eg_T_Gam = GaP_Eg0_Gam + GaP_alp_Gam * (1 - coth(GaP_Beta_Gam/T));
GaP_Eg_T_X = GaP_Eg0_X - GaP_alp_X * T^2 / (T + GaP_Beta_X);
GaP_Eg_T_L = GaP_Eg0_L - GaP_alp_L * T^2 / (T + GaP_Beta_L);

GaAs_Eg_T_Gam = GaAs_Eg0_Gam - GaAs_alp_Gam * T^2 / (T + GaAs_Beta_Gam);
GaAs_Eg_T_X = GaAs_Eg0_X - GaAs_alp_X * T;
GaAs_Eg_T_L = GaAs_Eg0_L - GaAs_alp_L * T^2 / (T + GaAs_Beta_L);


    GaAsP_Eg_Gam = (1-GaAsP_comp) * GaP_Eg_T_Gam + GaAsP_comp* GaAs_Eg_T_Gam ...
        - GaAsP_comp * (1-GaAsP_comp) * GaAsP_C_Gam
    GaAsP_Eg_X = (1-GaAsP_comp) .* GaP_Eg_T_X + GaAsP_comp .* GaAs_Eg_T_X ...
        - GaAsP_comp * (1-GaAsP_comp) * GaAsP_C_X
    GaAsP_Eg_L = (1-GaAsP_comp) .* GaP_Eg_T_L + GaAsP_comp .* GaAs_Eg_T_L ...
        - GaAsP_comp * (1-GaAsP_comp) * GaAsP_C_L
    GaAsP_Eg_compare=[GaAsP_Eg_Gam, GaAsP_Eg_X, GaAsP_Eg_L]
    GaAsP_Eg_T = min(GaAsP_Eg_compare);
   

Eg=GaAsP_Eg_T
Ec=Eg/2; %[eV]
Ev=-Ec; %[eV]
Ed=Ec;
Ea=Ec;
%% Effective Mass Calculation
Fcomp=F0+(F1-F0)*GaAsP_comp+Fbow*GaAsP_comp*(1-GaAsP_comp);
Epcomp=Ep0+(Ep1-Ep0)*GaAsP_comp+Epbow*GaAsP_comp*(1-GaAsP_comp);
dSOcomp=dSO0+(dSO1-dSO0)*GaAsP_comp+Epbow*GaAsP_comp*(1-GaAsP_comp);
meeff=1/((1+2*Fcomp)+Epcomp*(Eg+2*dSOcomp/3)/(Eg*(Eg+dSOcomp)));
me=meeff*m0;
mheff=mh0+(mh1-mh0)*GaAsP_comp;
mh=mheff*m0;
%% Density of States, Ni, and Fermi Level
Nc=2*(2*pi*meeff*m0*kb*T/(h^2))^(3/2)*m_to_cm;
Nv=2*(2*pi*mheff*m0*kb*T/(h^2))^(3/2)*m_to_cm;
ni=((Nc*Nv)^.5)*exp((-Eg)./(2*Vt));

Efn=fsolve(@(Efn) fermi_charge_neutrality(Efn,Ec,k,T,Ev,Nc,Nv,0,ga,Ea,N,gd,Ed),Ec-1e-8,optimset('Display','off')) %[eV]
%% Mobility
Tn=T/300;

beta=24*epsilon_0*erel*me*(k*T*q)^2/(hbar^2*q^2*N*cm_to_m);







mumaxp=150; %[cm2/Vs]
Nrefp=6.25e17; %[cm-3]
erel=11.1+1.8*GaAsP_comp;%relative dielectric permittivity
mu=mumaxp*(Tn.^(-2.3))./(1+N/(Nrefp*Tn^3.8)) %nGaAs base (minority holes in n type)
