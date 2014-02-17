%% Initialization
%=================================================
function [mu]=mobility(T,N,Cl,me,mh, epsilon_rel, debyeT, E1, K, pn)
%% Inputs
%Cl
%N doping
%me,mh  effective mass
%mh
%doping
%E1 deformation potential
%K electromechanical constant
%Eg
%E1=10.66; %[eV] (Deformation potential -- average of 21 experimental and calculated values from Adachi) 10.66
%K=.0617; % Electromechanical coupling constant -- Adachi) Ksq=3.7e-4;
%% Constants
%=================================================
k=8.6173324e-5; %[eV/K] (Wikipedia)
h=6.62606957e-34; %[J.s] (Wikipedia)
hbar=h/(2*pi); %[J.s] (obvs)
q=1.602176565e-19; %[Coul] (Wikipedia)
m=9.10938215e-31; %[kg] (Wikipedia)
Zd=1; %Absolute value of ion charge of acceptor
Za=1; %Absolute value of ion charge of donor
epsilon_0=8.854187817620e-12; %[F/m]
Ed_Si=0.00001; %[eV] (Pierret)0.0058 
Ea_Zn=0.00001; %[eV] (Pierret) 0.031
m_to_cm=1e-6; % Conversion from m^-3 to cm^-3
cm_to_m=1e6; % Conversion fromc cm^-3 to m^-3
cm2_to_m2=1e-4; % Conversion from cm^2 to m^2
m2_to_cm2=1e4; % Conversion from m^2 to cm^2
dyn_to_N=1e-5; % Conversion from dyne to Newton
gd=2; %Dengeneracy of electrons in conduction band
ga=4; %Degeneracy of heavy holes and light holes in valence band
Nni=0;
%% Determination of mobility
alphac=22.88*sqrt(me)*sqrt(300/debyeT)*1/(epsilon_0*epsilon_rel); % Coupling constant (Boer)
beta=24.*epsilon_0.*epsilon_rel.*me.*(k.*T.*q).^2./(hbar.^2.*q.^2.*N.*cm_to_m); %additional q on top is a conversion factor from eV to J

mudp_n=2.*sqrt(2.*pi).*hbar.^4.*Cl.*dyn_to_N.*cm2_to_m2.^-1.*q./(3.*(E1.*q).^2.*(me).^(5./2).*(k.*T.*q).^(3./2)).*m2_to_cm2;%deformation potential
mupe_n=16.*sqrt(2.*pi).*epsilon_0.*epsilon_rel.*hbar.^2./(3.*q.*K.^2.*me.^(3./2).*(k.*T.*q).^(1./2)).*m2_to_cm2;%piezoelectric

%mupo_n28=q.*hbar./(2.*me.*alphac.*k.*q.*debyeT).*exp(debyeT./T).*m2_to_cm2;%polar-optical phonons
mupo_n30=(3./2).*q.*hbar./(2.*me.*alphac.*k.*q.*debyeT).*sqrt(T./debyeT).*exp(debyeT./T).*m2_to_cm2;
%mupo_n31
%mupo_n32
%% Electron mobility, p-layer
if pn=1
muni_Pn=q.^3.*me./(80.*pi.*(Nni).*hbar.^3.*epsilon_0.*epsilon_rel.*cm_to_m).*m2_to_cm2;
muii_Pn=128.*sqrt(2.*pi).*(epsilon_0.*epsilon_rel).^2.*(k.*T.*q).^(3./2)./(N.*cm_to_m.*Za.^2.*q.^3.*(me).^(1./2).*(log(1+beta)-beta./(1+beta))).*m2_to_cm2;
mu_Pn=(mudp_n.^-1+mupe_n.^-1+muni_Pn.^-1+muii_Pn.^-1+mupo_n30.^-1).^-1;
mu=double(mu_Pn);
%Hole mobility, p-layer
else
%Electron mobility, n-layer
muni_Nn=q.^3.*me./(80.*pi.*(Nni).*hbar.^3.*epsilon_0.*epsilon_rel.*cm_to_m).*m2_to_cm2; %+Nd-Nd_ion(a
muii_Nn=128.*sqrt(2.*pi).*(epsilon_0.*epsilon_rel).^2.*(k.*T.*q).^(3./2)./(N.*cm_to_m.*Zd.^2.*q.^3.*(me).^(1./2).*(log(1+beta)-beta./(1+beta))).*m2_to_cm2;
mu_Nn=(mudp_n.^-1+mupe_n.^-1+muni_Nn.^-1+muii_Nn.^-1+(mupo_n30).^-1).^-1;
%Hole mobility, n-layer
end
%%
%===================
%Carrier Transport
%===================
%Holes
%vth_p=(3.*k.*T./mh).^(1/2); %[cm/s]
%taup=1./(vth_p.*sigmah.*Nt); %[s]

%electrons
%vth_n=(3.*k.*T./me).^(1/2); %[cm/s]
%taun=1/(vth_n.*sigmap.*Nt); %[s]
