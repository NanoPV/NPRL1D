function [Eg, ni, mu, erel]=nGaAsP_Eg(GaAsP_comp, N, T)
%In xAs
%GaP_0
%GaAs_1
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
%% Physical Parameters
c11_0=14.05E11; %[dyne/cm2]
c11_1=11.88E11;
c12_0=6.2E11;
c12_1=5.34E11;
c44_0=7.03E11;
c44_1=5.94E11;
K %electromechanical coupling constant

%% DeBye Temperature
M2=69.723;%Ga Atomic mass
M1=74.922;%As Atomic mass
M0=30.974;%P Atomic Mass
Mavg=(M2+M0+(M1-M0)*GaAsP_comp)/2


alphac

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
GaP_Eg0_Gam = 2.886;
GaP_alp_Gam = 0.1081;
GaP_Beta_Gam = 164;
GaP_Eg0_X = 2.35;
GaP_alp_X = 0.5771e-3;
GaP_Beta_X = 372;
GaP_Eg0_L = 2.72;
GaP_alp_L = 0.5771e-3;
GaP_Beta_L = 372;

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
erel=11.1+1.8*GaAsP_comp;%relative dielectric permittivity
%% Effective Mass Calculation
Fcomp=F0+(F1-F0)*GaAsP_comp+Fbow*GaAsP_comp*(1-GaAsP_comp);
Epcomp=Ep0+(Ep1-Ep0)*GaAsP_comp+Epbow*GaAsP_comp*(1-GaAsP_comp);
dSOcomp=dSO0+(dSO1-dSO0)*GaAsP_comp+Epbow*GaAsP_comp*(1-GaAsP_comp);
meeff=1/((1+2*Fcomp)+Epcomp*(Eg+2*dSOcomp/3)/(Eg*(Eg+dSOcomp)));
me=meeff*m0;
mheff=mh0+(mh1-mh0)*GaAsP_comp;
mh=mheff*m0;

%% Elasticity
c11=c11_0+(c11_1-c11_0)*GaAsP_comp;
c12=c12_0+(c12_1-c12_0)*GaAsP_comp;
c44=c44_0+(c44_1-c44_0)*GaAsP_comp;
Cl_100=c11; %[dyn/cm^2]
Cl_110=(1/2)*(c11+c12+c44); %[dyn/cm^2]
Cl_111=(1/5)*(3*c11+2*c12+4*c44); %[dyn/cm^2]
Cl=(Cl_100+Cl_110+Cl_111)/3; %[dyn/cm^2] 

%% Density of States, Ni, and Fermi Level
Nc=2*(2*pi*meeff*m0*kb*T/(h^2))^(3/2)*m_to_cm;
Nv=2*(2*pi*mheff*m0*kb*T/(h^2))^(3/2)*m_to_cm;
ni=((Nc*Nv)^.5)*exp((-Eg)./(2*Vt));
%% Debye tempetature
debyeT0=344;
debyeT1=445;

syms (debyeB debyeM)
debyeMatrix=solve(debyeB+debyeM*(M0+M2)/2==debyeT0, debyeB+debyeM*(69.723+30.974)/2==debyeT1);
debyeT=debyeMatrix.debyeB+debyeMatrix.debyeM*Mavg;



%% Mobility
mu=mobility(T,N,Cl,me,mh,erel,debyeT,E1,K,pn)

