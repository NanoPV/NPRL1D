function [Eg,ni]=AlAsSb_Eg(AlAsSb_comp,T)
%In xSB
%AlAs_0
%AlSb_1
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

AlSb_Eg0_Gam = 2.386; %[eV]
AlSb_alp_Gam = 0.42e-3; %[eV/K]
AlSb_Beta_Gam = 140; %[K]
AlSb_Eg0_X = 1.696; %[eV]
AlSb_alp_X = 0.39e-3; %[eV/K]
AlSb_Beta_X = 140; %[K]
AlSb_Eg0_L = 2.329; %[eV]
AlSb_alp_L = 0.58e-3; %[eV/K]
AlSb_Beta_L = 140; %[K]
AlSb_Lat_Const = 6.1355 + 2.60e-5 * (T-300); %[Ang]

AlAsSb_C_Gam = 0.8;
AlAsSb_C_X = 0.28;
AlAsSb_C_L = 0.28;
AlAsSb_Lat_Const = AlAs_Lat_Const + AlAsSb_comp * abs(AlAs_Lat_Const - AlSb_Lat_Const);

AlAs_Eg_T_Gam = AlAs_Eg0_Gam + AlAs_alp_Gam * (1 - coth(AlAs_Beta_Gam/T));
AlAs_Eg_T_X = AlAs_Eg0_X - AlAs_alp_X * T^2 / (T + AlAs_Beta_X);
AlAs_Eg_T_L = AlAs_Eg0_L - AlAs_alp_L * T^2 / (T + AlAs_Beta_L);

AlSb_Eg_T_Gam = AlSb_Eg0_Gam - AlSb_alp_Gam * T^2 / (T + AlSb_Beta_Gam);
AlSb_Eg_T_X = AlSb_Eg0_X - AlSb_alp_X * T;
AlSb_Eg_T_L = AlSb_Eg0_L - AlSb_alp_L * T^2 / (T + AlSb_Beta_L);


    AlAsSb_Eg_Gam = (1-AlAsSb_comp) * AlAs_Eg_T_Gam + AlAsSb_comp* AlSb_Eg_T_Gam ...
        - AlAsSb_comp * (1-AlAsSb_comp) * AlAsSb_C_Gam
    AlAsSb_Eg_X = (1-AlAsSb_comp) .* AlAs_Eg_T_X + AlAsSb_comp .* AlSb_Eg_T_X ...
        - AlAsSb_comp * (1-AlAsSb_comp) * AlAsSb_C_X
    AlAsSb_Eg_L = (1-AlAsSb_comp) .* AlAs_Eg_T_L + AlAsSb_comp .* AlSb_Eg_T_L ...
        - AlAsSb_comp * (1-AlAsSb_comp) * AlAsSb_C_L
    AlAsSb_Eg_compare=[AlAsSb_Eg_Gam, AlAsSb_Eg_X, AlAsSb_Eg_L];
    AlAsSb_Eg_T = min(AlAsSb_Eg_compare);
    
 Eg=AlAsSb_Eg_T;
