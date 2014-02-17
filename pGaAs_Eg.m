function [Eg_GaAs, ni_GaAs, mu, erel]=pGaAs_Eg(Na, Nd, T)

kb=1.381e-23;
q=1.602000000000000e-19;
erel=13.1;%erel_GaAs;




Vt=kb*T/q;
Eg_GaAs=1.519-5.405*10^-4*T^2/(204+T); %in eV, for intrinsic GaAs  
Nc_GaAs=8.63e13*T^(3/2)*(1-1.93e-4*T-4.19e-8*T^2);
Nv_GaAs=1.83e15*T^(3/2);
dEg_GaAs=2e-11*Ne_a.^(.5); %in eV, p-type gap narrowing from Marti
dEgn_GaAs=.125/erel*abs(1-Nb_d./Ncr).^(1/3); %eV (intrinsic (n-type?) gap narrowing from Marti)
ni_GaAs=((Nc_GaAs*Nv_GaAs)^.5)*exp((-Eg_GaAs-dEgn_GaAs)./(2*Vt)); %ni in base

c=1/3;
mumaxn=-6933*c^3+3581*c^2+76*c+9157;
muminn=-2201*c^3+5362*c^2-5092*c+1942;
Nrefn=-4.36e16*c^3+1.43e17*c^2-1.80e17*c+8.49e16;
chi=.506*c^3-.262*c^2+.132*c+.365;
Tn=T/300;

mumaxp=150; %[cm2/Vs]
Nrefp=6.25e17; %[cm-3]

mu=muminn*Tn.^(-.57)+((mumaxn-muminn)*Tn^(-2.3))./(1+(Ma./(Nrefn*Tn^(-2.4))).^(chi*Tn^(-.146))); %pGaAs emitter (minority e- in p type)
