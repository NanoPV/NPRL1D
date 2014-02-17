function [Eg, ni, mu, meeff,mheff, erel]=InAlAs_Eg(comp, N, T)
%In xAs
%AlAs_0
%InAs_1
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
%c11_0=14.05E11; %[dyne/cm2]
%c11_1=11.88E11;
%c12_0=6.2E11;
%c12_1=5.34E11;
%c44_0=7.03E11;
%c44_1=5.94E11;
%Cl_100=c11; %[dyn/cm^2]
%Cl_110=(1/2)*(c11+c12+c44); %[dyn/cm^2]
%Cl_111=(1/5)*(3*c11+2*c12+4*c44); %[dyn/cm^2]
%Cl=(Cl_100+Cl_110+Cl_111)/3; %[dyn/cm^2] 


%K %electromechanical coupling constant

%% DeBye Temperature
M2=69.723;%Ga Atomic mass
M1=74.922;%As Atomic mass
M0=30.974;%P Atomic Mass
Mavg=(M2+M0+(M1-M0)*comp)/2;


%alphac

%% Bowing Parameters for effective mass
F0=-0.48;
F1=-1.94;
Fbow=-4.44;

Ep0=21.1;
Ep1=28.8;
Epbow=-4.81;

dSO0=0.28;
dSO1=0.341;
dSObow=0.15;

mh0=0.8;
mh1=0.41;
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
AlAs_Lat_Const = 5.6611 + 2.9e-5 * (T-300); %[Ang]

InAs_Eg0_Gam = 0.417; %[eV]
InAs_alp_Gam = 0.276e-3; %[eV/K]
InAs_Beta_Gam = 93; %[K]
InAs_Eg0_X = 1.433; %[eV]
InAs_alp_X = 0.276e-3; %[eV/K]
InAs_Beta_X = 93; %[K]
InAs_Eg0_L = 1.133; %[eV]
InAs_alp_L = 0.276e-3; %[eV/K]
InAs_Beta_L = 93; %[K]
InAs_Lat_Const = 6.0583 + 2.74e-5 * (T-300); %[Ang]

AlInAs_C_Gam = 0.70;
AlInAs_C_X = 0;
AlInAs_C_L = 0;
AlInAs_Lat_Const = AlAs_Lat_Const + comp * abs(AlAs_Lat_Const - InAs_Lat_Const);


AlAs_Eg_T_Gam = AlAs_Eg0_Gam + AlAs_alp_Gam * (1 - coth(AlAs_Beta_Gam/T));
AlAs_Eg_T_X = AlAs_Eg0_X - AlAs_alp_X * T^2 / (T + AlAs_Beta_X);
AlAs_Eg_T_L = AlAs_Eg0_L - AlAs_alp_L * T^2 / (T + AlAs_Beta_L);

InAs_Eg_T_Gam = InAs_Eg0_Gam - InAs_alp_Gam * T^2 / (T + InAs_Beta_Gam);
InAs_Eg_T_X = InAs_Eg0_X - InAs_alp_X * T;
InAs_Eg_T_L = InAs_Eg0_L - InAs_alp_L * T^2 / (T + InAs_Beta_L);


    AlInAs_Eg_Gam = (1-comp) * AlAs_Eg_T_Gam + comp* InAs_Eg_T_Gam ...
        - comp * (1-comp) * AlInAs_C_Gam
    AlInAs_Eg_X = (1-comp) .* AlAs_Eg_T_X + comp .* InAs_Eg_T_X ...
        - comp * (1-comp) * AlInAs_C_X
    AlInAs_Eg_L = (1-comp) .* AlAs_Eg_T_L + comp .* InAs_Eg_T_L ...
        - comp * (1-comp) * AlInAs_C_L
    AlInAs_Eg_compare=[AlInAs_Eg_Gam, AlInAs_Eg_X, AlInAs_Eg_L]
    AlInAs_Eg_T = min(AlInAs_Eg_compare);
   

Eg=AlInAs_Eg_T
Ec=Eg/2; %[eV]
Ev=-Ec; %[eV]
Ed=Ec;
Ea=Ec;
erel=11.1+1.8*comp;%relative dielectric permittivity
%% Effective Mass Calculation
Fcomp=F0+(F1-F0)*comp+Fbow*comp*(1-comp);
Epcomp=Ep0+(Ep1-Ep0)*comp+Epbow*comp*(1-comp);
dSOcomp=dSO0+(dSO1-dSO0)*comp+Epbow*comp*(1-comp);
meeff=1/((1+2*Fcomp)+Epcomp*(Eg+2*dSOcomp/3)/(Eg*(Eg+dSOcomp)));
meeff;
me=meeff*m0;
mheff=mh0+(mh1-mh0)*comp;
mh=mheff*m0;
mheff;
%% Elasticity
%c11=c11_0+(c11_1-c11_0)*comp;
%c12=c12_0+(c12_1-c12_0)*comp;
%c44=c44_0+(c44_1-c44_0)*comp;
%Cl_100=c11; %[dyn/cm^2]
%Cl_110=(1/2)*(c11+c12+c44); %[dyn/cm^2]
%Cl_111=(1/5)*(3*c11+2*c12+4*c44); %[dyn/cm^2]
%Cl=(Cl_100+Cl_110+Cl_111)/3; %[dyn/cm^2] 


%% Density of States, Ni, and Fermi Level
Nc=2*(2*pi*meeff*m0*kb*T/(h^2))^(3/2)*m_to_cm;
Nv=2*(2*pi*mheff*m0*kb*T/(h^2))^(3/2)*m_to_cm;
ni=((Nc*Nv)^.5)*exp((-Eg)./(2*Vt))
%% Debye tempetature
debyeT0=344;
debyeT1=445;

syms debyeB debyeM
debyeMatrix=solve(debyeB+debyeM*(M0+M2)/2==debyeT0, debyeB+debyeM*(M1+M2)/2==debyeT1);
debyeT=debyeMatrix.debyeB+debyeMatrix.debyeM*Mavg;



%% Mobility

mu=1;






