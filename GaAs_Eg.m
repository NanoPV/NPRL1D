function [Eg_pGaAs, ni_GaAs, mue, mui, mub, musub]=GaAs_Eg(Ne_a, Nb_d, T)

kb=1.381e-23;
q=1.602000000000000e-19;
erel_e=13.1;%erel_GaAs;
erel_i=13.1;%erel_GaAs; 
erel_b=13.1;%erel_GaAs; 
Ncr=1.6e24*(.066/(1.4*erel_b))^3;
Maw=Naw; %maj carrier conc in window
Ma=Na; %maj carrier conc in emitter
Mi=Ni;
Md=Nd; %maj carrier conc in base
Msub=Nsub;


Vt=kb*T/q;
Eg_GaAs=1.519-5.405*10^-4*T^2/(204+T); %in eV, for intrinsic GaAs  
Nc_GaAs=8.63e13*T^(3/2)*(1-1.93e-4*T-4.19e-8*T^2);
Nv_GaAs=1.83e15*T^(3/2);
dEg_GaAs=2e-11*Ne_a.^(.5); %in eV, p-type gap narrowing from Marti
dEgn_GaAs=.125/erel_e*abs(1-Nb_d./Ncr).^(1/3); %eV (intrinsic (n-type?) gap narrowing from Marti)
ni_GaAs=((Nc_GaAs*Nv_GaAs)^.5)*exp((-Eg_GaAs-dEgn_GaAs)./(2*Vt)); %ni in base

c=1/3;
mumaxn=-6933*c^3+3581*c^2+76*c+9157;
muminn=-2201*c^3+5362*c^2-5092*c+1942;
Nrefn=-4.36e16*c^3+1.43e17*c^2-1.80e17*c+8.49e16;
chi=.506*c^3-.262*c^2+.132*c+.365;
Tn=T/300;

mumaxp=150; %[cm2/Vs]
Nrefp=6.25e17; %[cm-3]

muw=muminn*Tn^(-.57)+((mumaxn-muminn)*Tn^(-2.3))/(1+(Maw/(Nrefn*Tn^(-2.4)))^(chi*Tn^(-.146))); %InGaP window  (minority electrons in ptype)
mue=muminn*Tn.^(-.57)+((mumaxn-muminn)*Tn^(-2.3))./(1+(Ma./(Nrefn*Tn^(-2.4))).^(chi*Tn^(-.146))); %pGaAs emitter (minority e- in p type)
mui=muminn*Tn.^(-.57)+((mumaxn-muminn)*Tn^(-2.3))./(1+(Mi./(Nrefn*Tn^(-2.4))).^(chi*Tn^(-.146))); %pGaAs i-region
mub=mumaxp*(Tn.^(-2.3))./(1+Md./(Nrefp*Tn^3.8)) %nGaAs base (minority holes in n type)
musub=1.5E4;%mumaxp*(Tn^(-2.3))/(1+Msub/(Nrefp*Tn^3.8));