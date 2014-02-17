function [Ef]=fermisolve(Eg,N, me, mh, pn, T)
%% Inputs
% Eg Bandgap[eV]
% me electron effective mass [kg]
% mh hole effective mass [kg]s
% pn 1 if p-type 2 if n-type
%% Outputs
%Ef Fermi Level: distance from respective band [eV]
%% Constants
k=8.6173324e-5; %[eV/K] (Wikipedia)
h=6.62606957e-34; %[J.s] (Wikipedia)
hbar=h/(2*pi); %[J.s] (obvs)
q=1.602176565e-19; %[Coul] (Wikipedia)
m=9.10938215e-31; %[kg] (Wikipedia)
Zd=1; %Absolute value of ion charge of acceptor
Za=1; %Absolute value of ion charge of donor
epsilon_0=8.854187817620e-12; %[F/m]; %(Boer)
Ed_Si=0.00001; %[eV] (Pierret)0.0058 
Ea_Zn=0.00001; %[eV] (Pierret) 0.031
m_to_cm=1e-6; % Conversion from m^-3 to cm^-3
cm_to_m=1e6; % Conversion fromc cm^-3 to m^-3
cm2_to_m2=1e-4; % Conversion from cm^2 to m^2
m2_to_cm2=1e4; % Conversion from m^2 to cm^2
dyn_to_N=1e-5; % Conversion from dyne to Newton
gd=2; %Dengeneracy of electrons in conduction band
ga=4; %Degeneracy of heavy holes and light holes in valence band
m0=9.11E-31; %kg
%% Determination of temperature dependent parameters 
Nc=2*((2*pi*me*m0*k*T*q)/(h^2)).^(3/2)*m_to_cm; %[cm^-3]
Nv=2*((2*pi*mh*m0*k*T*q)/(h^2)).^(3/2)*m_to_cm; %[cm^-3]
Ec=Eg/2; %[eV]
Ev=-Ec; %[eV]
Edb=Ec-3*k*T; %[eV]
Eab=Ev+3*k*T; %[eV]
Ed=Ec-Ed_Si; %[eV]
Ea=Ev+Ea_Zn; %[eV]
Ei=(Ec+Ev)/2+(k*T)/2*log(Nv/(Nc)); %[eV]
%% Determination of Ef
%Emitter
if pn==1
Efp=fsolve(@(Efp) fermi_charge_neutrality(Efp,Ec,k,T,Ev,Nc,Nv,N,ga,Ea,0,gd,Ed),Ev+1e-8,optimset('Display','off')); %[eV]
etacp=(Efp-Ec)/(k*T);
etavp=(Ev-Efp)/(k*T);
fermiIntcp=@(x) sqrt(x)./(1+exp(x-etacp));
fermiIntvp=@(x) sqrt(x)./(1+exp(x-etavp));
fermicp=1/(gamma(3/2))*quadgk(fermiIntcp,1E-31,20);
fermivp=1/(gamma(3/2))*quadgk(fermiIntvp,0,20);
Pn=Nc*fermicp; %[cm^-3] minotiry electrons
Pp=Nv*fermivp; %[cm^-3] majority holes
Na_ion =N/(1+ga*exp((Ea -Efp )/(k*T ))); %[cm^-3]
Ef=Efp;
%Base
else
Efn=fsolve(@(Efn) fermi_charge_neutrality(Efn,Ec,k,T,Ev,Nc,Nv,0,ga,Ea,N,gd,Ed),Ec-1e-8,optimset('Display','off')) %[eV]%Determination of carrier densities
etacn=(Efn -Ec )/(k*T );
etavn=(Ev -Efn )/(k*T );
fermiIntcn=@(x) sqrt(x)./(1+exp(x-etacn));
fermiIntvn=@(x) sqrt(x)./(1+exp(x-etavn));
fermicn=1/(gamma(3/2))*quadgk(fermiIntcn,0,10);
fermivn=1/(gamma(3/2))*quadgk(fermiIntvn,0,inf);
Nn=Nc*fermicn; %[cm^-3] majority electrons
Np=Nv*fermivn; %[cm^-3] minority holes
Nd_ion =N/(1+gd*exp((Efn-Ed)/(k*T))); %[cm^-3]
Ef=Efn;
end