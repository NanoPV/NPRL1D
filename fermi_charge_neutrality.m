function F = fermi_charge_neutrality(Ef,Ec,k,T,Ev,Nc,Nv,Na,ga,Ea,Nd,gd,Ed)
%Charge neutrality invoking Fermi-Dirac statistics
%Equations from Pierret, Advanced Semiconductor Fundamentals, 2nd Ed, 2003.
etac=(Ef-Ec)./(k.*T);
etav=(Ev-Ef)./(k.*T);
fermiIntc=@(x) sqrt(x)./(1+exp(x-etac));
fermiIntv=@(x) sqrt(x)./(1+exp(x-etav));
fermic=1/(gamma(3/2)).*quadgk(fermiIntc,0,inf);
fermiv=1/(gamma(3/2)).*quadgk(fermiIntv,0,inf);
n=Nc*fermic;
p=Nv*fermiv;
Na_ion=Na./(1+ga*exp((Ea-Ef)./(k*T)));
Nd_ion=Nd./(1+gd*exp((Ef-Ed)./(k*T)));
F=p-n+Nd_ion-Na_ion;
end