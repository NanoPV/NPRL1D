%% Initialization
%=================================================
function [mup, mun]=mobility(Cl,me,mh, epsilon_rel, debyeT, alphac, Eg
%% Required parameters
%Cl 
%me,mh  effective mass
%mh
%doping
%E1 deformation potential
%K electromechanical constant
%Eg

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
epsilon_rel=12.9; %(Boer)
Ed_Si=0.00001; %[eV] (Pierret)0.0058 
Ea_Zn=0.00001; %[eV] (Pierret) 0.031
m_to_cm=1e-6; % Conversion from m^-3 to cm^-3
cm_to_m=1e6; % Conversion fromc cm^-3 to m^-3
cm2_to_m2=1e-4; % Conversion from cm^2 to m^2
m2_to_cm2=1e4; % Conversion from m^2 to cm^2
dyn_to_N=1e-5; % Conversion from dyne to Newton
gd=2; %Dengeneracy of electrons in conduction band
ga=4; %Degeneracy of heavy holes and light holes in valence band
%% Material Variables
%=================================================

E1=10.66; %[eV] (Deformation potential -- average of 21 experimental and calculated values from Adachi) 10.66
K=.0617; % Electromechanical coupling constant -- Adachi) Ksq=3.7e-4;

%% Determination of temperature dependent parameters 
%=================================================    
Nc=2*((2*pi*me*k*T*q)/(h^2)).^(3/2)*m_to_cm; %[cm^-3]
Nv=2*((2*pi*mh*k*T*q)/(h^2)).^(3/2)*m_to_cm; %[cm^-3]
Ec=Eg/2; %[eV]
Ev=-Ec; %[eV]
Edb=Ec-3*k*T; %[eV]
Eab=Ev+3*k*T; %[eV]
Ed=Ec-Ed_Si; %[eV]
Ea=Ev+Ea_Zn; %[eV]
Ei=(Ec+Ev)/2+(k*T)/2*log(Nv/(Nc)); %[eV]
ni=sqrt(Nc.*Nv).*exp(-(Eg)./(2*k*T)); %[cm^-3]

%% Determination of Ef


%Base

%=================================================
%Emitter
Efp=fsolve(@(Efp) fermi_charge_neutrality(Efp,Ec,k,T,Ev,Nc,Nv,Na,ga,Ea,0,gd,Ed),Ev+1e-8,optimset('Display','off')); %[eV]
etacp=(Efp-Ec)/(k*T);
etavp=(Ev-Efp)/(k*T);
fermiIntcp=@(x) sqrt(x)/(1+exp(x-etacp));
fermiIntvp=@(x) sqrt(x)/(1+exp(x-etavp));
fermicp=1/(gamma(3/2))*quadgk(fermiIntcp,0,inf);
fermivp=1/(gamma(3/2))*quadgk(fermiIntvp,0,inf);
Pn=Nc*fermicp; %[cm^-3] minotiry electrons
Pp=Nv*fermivp; %[cm^-3] majority holes
Na_ion =Na/(1+ga*exp((Ea -Efp )/(k*T ))); %[cm^-3]
%Base
Efn=fsolve(@(Efn) fermi_charge_neutrality(Efn,Ec,k,T,Ev,Nc,Nv,0,ga,Ea,N,gd,Ed),Ec-1e-8,optimset('Display','off')) %[eV]%Determination of carrier densities
etacn=(Efn -Ec )/(k*T );
etavn=(Ev -Efn )/(k*T );
fermiIntcn=@(x) sqrt(x)./(1+exp(x-etacn));
fermiIntvn=@(x) sqrt(x)./(1+exp(x-etavn));
fermicn=1/(gamma(3/2))*quadgk(fermiIntcn,0,inf);
fermivn=1/(gamma(3/2))*quadgk(fermiIntvn,0,inf);
Nn=Nc*fermicn; %[cm^-3] majority electrons
Np=Nv*fermivn; %[cm^-3] minority holes
Nd_ion =Nd/(1+gd*exp((Efn-Ed)/(k*T))); %[cm^-3]
%%
eta=q*hbar./(2*k.*T*q).*sqrt(Na.*cm_to_m./(epsilon_0*epsilon_rel*mh));
n=eta.*coth(eta);
plot(T,n);
%%
Eckt=Ec-(3*k.*T);
Evkt=Ev+(3*k.*T);
plot(T,Ec,T,Ev,T,Efn,T,Efp,T,Eckt,T,Evkt)

%% Determination of mobility
%=================================================

alphac=22.88*sqrt(me)*sqrt(300/debyeT)*1/(epsilon_0*epsilon_gaas); % Coupling constant (Boer)
beta=24.*epsilon_0.*epsilon_rel.*me.*(k.*T.*q).^2./(hbar.^2.*q.^2.*Nn.*cm_to_m); %additional q on top is a conversion factor from eV to J

mudp_n=2.*sqrt(2.*pi).*hbar.^4.*Cl.*dyn_to_N.*cm2_to_m2.^-1.*q./(3.*(E1.*q).^2.*(me).^(5./2).*(k.*T.*q).^(3./2)).*m2_to_cm2;
mupe_n=16.*sqrt(2.*pi).*epsilon_0.*epsilon_rel.*hbar.^2./(3.*q.*K.^2.*me.^(3./2).*(k.*T.*q).^(1./2)).*m2_to_cm2;

mupo_n28=q.*hbar./(2.*me.*alphac.*k.*q.*debyeT).*exp(debyeT./T).*m2_to_cm2;
mupo_n30=(3./2).*q.*hbar./(2.*me.*alphac.*k.*q.*debyeT).*sqrt(T./debyeT).*exp(debyeT./T).*m2_to_cm2;
%mupo_n31
%mupo_n32

%% Electron mobility, p-emitter
muni_Pn=q.^3.*me./(80.*pi.*(Nni).*hbar.^3.*epsilon_0.*epsilon_rel.*cm_to_m).*m2_to_cm2;
muii_Pn=128.*sqrt(2.*pi).*(epsilon_0.*epsilon_rel).^2.*(k.*T.*q).^(3./2)./(Na_ion.*cm_to_m.*Za.^2.*q.^3.*(me).^(1./2).*(log(1+beta)-beta./(1+beta))).*m2_to_cm2;
mu_Pn=(mudp_n.^-1+mupe_n.^-1+muni_Pn.^-1+muii_Pn.^-1+mupo_n30.^-1).^-1;
%Hole mobility, p-emitter
%Electron mobility, n-base
muni_Nn=q.^3.*me./(80.*pi.*(Nni).*hbar.^3.*epsilon_0.*epsilon_rel.*cm_to_m).*m2_to_cm2; %+Nd-Nd_ion(a
muii_Nn=128.*sqrt(2.*pi).*(epsilon_0.*epsilon_rel).^2.*(k.*T.*q).^(3./2)./(Nd_ion.*cm_to_m.*Zd.^2.*q.^3.*(me).^(1./2).*(log(1+beta)-beta./(1+beta))).*m2_to_cm2;
mu_Nn=(mudp_n.^-1+mupe_n.^-1+muni_Nn.^-1+muii_Nn.^-1+(mupo_n30).^-1).^-1;
%Hole mobility, n-base
%waitbar(a/length(T),progress)
%end
%close(progress)

loglog(T,mudp_n,T,mupe_n,T,muni_Pn,T,muii_Pn,T,mu_Pn,T,mupo_n30)
axis([1 1e4 1e3 1e8]);
legend('mudp','mupe','muni','muii','mu','mupo_n30');

%%
%===================
%Carrier Transport
%===================
%Holes
%vth_p=(3.*k.*T./mh).^(1/2); %[cm/s]
%taup=1./(vth_p.*sigmah.*Nt); %[s]
%Dp=mup.*k.*T./q; %[cm^2/s]
%Lp=(Dp.*taup).^(1/2); %[cm]

%electrons
%vth_n=(3.*k.*T./me).^(1/2); %[cm/s]
%taun=1/(vth_n.*sigmap.*Nt); %[s]
%Dn=mun.*k.*T./q; %[cm^2/s]
%Ln=(Dn.*taun).^(1/2); %[cm]

%Jsat=q*peq*Lp/taup;
%% Electrostatics




