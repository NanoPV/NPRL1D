%NPRL Solar Cell Model
%Michael Harris, Christopher Kerestes, Zachary Bittner

%function [ Jsc, Voc, FFs, Efficiency ] = QEModelingGaP(X)

clear all
clc
close all
pause(1)

load('~/Google Drive/NPRL1D/QEworkspaceGaP.mat'); %absorption, spectra, and reflectance data
%WavelengthEQEIQERSRdata = xlsread('C:\Users\Zac Bittner\Desktop\Kerestes QE\GaAs60XR7.xls'); %file location of spectral response/EQE data
%Wave = WavelengthEQEIQERSRdata(:,1); EQE=WavelengthEQEIQERSRdata(:,3); SR=WavelengthEQEIQERSRdata(:,2);

%% incident solar spectra W/cm^2/nm for am0, am1.5g, amd1.5d 280-1100nm
spect=3; %spectrum to use [1,2,3] = [AM0,AM1.5G,AM1.5D]
F=Fsun(1:end,spect);
Psun=Fsun(1:end,spect); %incident solar power W/cm^2 for am0, am1.5g, amd1.5d

%% Reflection
ref=.05;%Ref;
%ref=0;

%% Variables
T=300 ; %[Kelvin]
Vt=kb*T/q;

%dielectric constants erel_GaAs=13.1, erel_InGaP=11.8 GaP 11.1

% K=13.1; %gaas dielectric constant
%Layer Thicknesses
t_w = .15e-5; %cm
t_e = .150e-4 %cm
t_i = .109E-4; %cm 8.73
%i-region thickness [cm]
%0X: 100e-7; 5X: 133.25e-7; 10X: 205e-7; 20X: 335e-7; 40X: 604e-7;
%60X:873e-7; 100X: 1411e-7;
t_b = 2e-4; %cm
t_sub = 350e-4; %cm


%doping (pwindow/pemitter/nbase) [cm-3]
Ncr=1.6e24*(.066/(1.4*erel_b))^3;%????
Naw=2e18;
Ne_a = 1e18; %changed Na to Ne_a
% Ni_a = ;%doping in i-region not included in Harris' program
Nb_d = 2e17; %1e17 changed Nd to Nb_d
Nsub=1e18;

Sw=1E9; %1e8; %front surface InGaP/air interface  %%%% VARIABLE!!!!!!!!!!!!!!!
Se=1E5; %emitter/window interface recombination  %%%% VARIABLE!!!!!!!!!!!!!!!
Sb=10000; %base back surface recombination  %%%% VARIABLE!!!!!!!!!!!!!!!
dataS=[Sw,Se,Sb];

%minority carrier lifetimes
tauw = 2e-12;%s%2E-13
taue = 2.5E-10;%1.4E-11%3.2E-11
taub = 100e-9; %s3.5E-11%3.5E-11
tauip = 1.1E-15%1.1E-15 %1.3E-15%1.2E-15;
tauie = 2.6E-11 %2.0E-11%2.1E-11;
taui=tauip;
% tausub1

%% Mobility and Diffusivity
%mobilities [cm2/Vs]   from Marti paper
%muw=1000; %InGaP window  (minority electrons in ptype)
%mue=250; %pGaAs emitter (minority e- in p type)
%mub=1000; %nGaAs base (minority holes in n type)
%musub=mub
%equations from Marti/Algora1997 
%Electrons.......  as min carriers c=(1+cp/3-cp)=1/3 with no compensation , M=maj carrier concentration
c=1/3;
mumaxn=250;%-6933*c^3+3581*c^2+76*c+9157;
muminn=75;%-2201*c^3+5362*c^2-5092*c+1942;
Nrefn=-4.36e16*c^3+1.43e17*c^2-1.80e17*c+8.49e16;
chi=.506*c^3-.262*c^2+.132*c+.365;
Tn=T/300;
% 
% %Aux equation
% mun=mumin*Tn^-.57+((mumax-mumin)*Tn^-2.3)/(1+(M/(Nref*Tn^-2.4))^(chi*Tn-.146))
%Holes........
mumaxp=150; %[cm2/Vs]
Nrefp=6.25e17; %[cm-3]
%  %Hole aux equation
%mup=mumaxp*Tn^-2.3/(1+M/(Nrefp*Tn^3.8))

muw=60;%muminn*Tn^(-.57)+((mumaxn-muminn)*Tn^(-2.3))/(1+(Maw/(Nrefn*Tn^(-2.4)))^(chi*Tn^(-.146))); %InGaP window  (minority electrons in ptype)
mue=130;%ioffe hole in InGaAs %muminn*Tn.^(-.57)+((mumaxn-muminn)*Tn^(-2.3))./(1+(Ma./(Nrefn*Tn^(-2.4))).^(chi*Tn^(-.146))) %pGaAs emitter (minority e- in p type)
mui=1.5E4;%muminn*Tn.^(-.57)+((mumaxn-muminn)*Tn^(-2.3))./(1+(Mi./(Nrefn*Tn^(-2.4))).^(chi*Tn^(-.146))); %pGaAs i-region
mub=130;%mumaxp*(Tn.^(-2.3))./(1+Md./(Nrefp*Tn^3.8)) %nGaAs base (minority holes in n type)
musub=1.5E4;%mumaxp*(Tn^(-2.3))/(1+Msub/(Nrefp*Tn^3.8));

%Diffusivity
Dw=Vt*muw;
De=Vt*mue
Di=Vt.*mui;
Db=Vt*mub

Lw=sqrt(Dw.*tauw)%cm
Le=sqrt(De.*taue)%cm
Li=2e-4; %cm 
Lb=sqrt(Db.*taub)%cm
dataL=[Lw,Le,Lb];


%% Temp and Concentration dependent Eg and ni options

%[Eg_InGaP, ni_InGaP, mu]=pInGaP_Eg(comp, Ne_a, Nb_d, T)
%[Eg_InGaP, ni_InGaP, mu]=nInGaP_Eg(comp, Ne_a, Nb_d, T)

%[Eg_GaAs, ni_GaAs, mue]=pGaAs_Eg(Ne_a, Nb_d, T);
%[Eg_GaAs, ni_GaAs, mui]=pGaAs_Eg(Ne_a, Nb_d, T);
%[Eg_GaAs, ni_GaAs, mub]=nGaAs_Eg(Ne_a, Nb_d, T);
xIn_InGaP=0.52;%LM to GaAs x=0.52

[Eg_InGaP, ni_InGaP, mue, erel_e]=pInGaP_Eg(xIn_InGaP, Ne_a, T)
[Eg_InGaP, ni_InGaP, mue, erel_i]=pInGaP_Eg(xIn_InGaP, ni_InGaP, T)
[Eg_InGaP, ni_InGaP, mub, erel_b]=pInGaP_Eg(xIn_InGaP, Nb_d, T)




GeComp=0.47;
Eg_InGaAs=0.42 + 0.625*GeComp -(5.8/(T+300)-4.19/(T+271))*1E-4*T^2*GeComp - 4.19E-4*T^2/(T+271)+0.475*GeComp^2;
Nc_InGaAs=4.82E15*(.023+.037*GeComp+.003*GeComp^2)^1.5*T^1.5;
Nv_InGaAs=4.82E15*(.41-.1*GeComp)^1.5*T^1.5;
ni_InGaAs=(Nc_InGaAs*Nv_InGaAs)^.5*exp(-Eg_InGaAs/(2*8.6173324E-5*T));%6.3E11

Eg_GaP=2.34-0.00061*(T^2)/(T+460);
Nc_GaP=3400000000000000*T^(3/2);
Nv_GaP=3600000000000000*T^(3/2);
ni_GaP=(Nc_GaP*Nv_GaP)^0.5*exp(-Eg_GaP/(0.00008617*T*2));

ni_w = 0.2;
ni_e = ni_InGaP; %Ga(x)In(1-x)P, x=51.5% Eg=1.85
ni_b = ni_InGaP;

%% built in junction voltage and depletion width
Vbi=Vt*log(Nb_d.*Ne_a./ni_e.^2);
Vbin=Vt*log(Nb_d/ni_b);
Vbip=Vt*log(Ne_a/ni_e);
Vbisep=Vt*( log(Nb_d./ni_b)+log(Ne_a./ni_e) );
Vbi=Vbisep;

% Wdi=(2*erel_i*e0/q*(Ne_a+Nb_d)./(Ne_a.*Nb_d).*(Vbi-.95)).^.5;
% WdiBart=(t_i.^2 + (2.*erel_i.*e0./q) .* (Ne_a+Nb_d)./(Ne_a.*Nb_d) .* (Vbi-0) ).^.5 - t_i; %from Bart
% % Wdi=WdiBart;
% Vpt = t_i.^2 * 1e15* q ./ (2*erel_i*e0);

wi_e=1./Ne_a.*sqrt(2*erel_e*e0*Vbi./(q*(1./Ne_a+1./Nb_d)));%changed wp to wi_e
wi_b=1./Nb_d.*sqrt(2*erel_b*e0*Vbi./(q*(1./Ne_a+1./Nb_d)));%changed wn to wi_b
Wd=wi_e+wi_b;


%% Absorption Coefficients
a_w=a_GaPAbsStaffan; %absorption in windowchanged Aw to a_w
a_e=a_In52Ga48P; %absorption in emitterchanged An to a_e
a_i=a_In52Ga48P; %absorption in i-regionGaAs changed Ai to a_i
a_b=a_In52Ga48P; %absorption in baseGaAs changed Ap to a_b

%subscript w=window (a in Hovel) e=emitter (g in Hovel), b=base (p Hovel)
%absorption coefficients taken from Aguinaldo data p=pGaAs, n=nGaAs, and
%measured GaAs data from Michael Slocum

%% Hovel Model 

%Constants for QE fitting in window layer
Q=a_w*Lw; %changed B to a_w
Y=-t_w/Lw; %changed D to t_w
Z=-Sw*tauw/Lw; %chaged tw to tauw

% Constants for QE fitting in emitter layer
S=a_e.*Le;
T2=-t_e/Le; %changed d to t_e
U=-Se*taue/Le; %changed te to taue
V = Se*Le/De;

% Constants for QE fitting in base layer
E=a_b.*Lb;
Fa=t_b/Lb; %changed Tb to t_b
G=Sb.*taub./Lb; %changed tb to taub
J = Sb.*Lb./Db;


%average absorption distance in i-region
%xavg=zeros(941,1);
%xavg=(-1./a_e.*log((1-exp(-a_e.*(t_i+Wd))./(a_e.*(t_i+Wd)))));
t_trans=(Wd+t_i)/1.2E7
%from Nelson:
pass=1
QE_w=(1-ref).*exp(-2.15E-4.*pass.*a_GaPAbsStaffan).*a_w.*Lw./(a_w.^2.*Lw.^2-1).*( (a_w.*Lw+Sw.*tauw./Lw.*(1-exp(-a_w.*t_w).*cosh(t_w./Lw))-exp(-a_w.*t_w).*sinh(t_w./Lw))./(Sw.*tauw./Lw.*sinh(t_w./Lw)+cosh(t_w./Lw))-a_w.*Lw.*exp(-a_w.*t_w) );
QE_e = (1-ref) .*exp(-2.15E-4.*pass.*a_GaPAbsStaffan).* exp(-a_w.*t_w) .* ( S./(S.^2 - 1)) .*  (  ((V+a_e.*Le) - exp(-a_e.*(t_e-wi_e)).*(V.*cosh((t_e-wi_e)./Le)+sinh((t_e-wi_e)./Le) ) ) / ( V.*sinh((t_e-wi_e)./Le)+cosh((t_e-wi_e)./Le) ) - (a_e.*Le.*exp(-a_e.*(t_e-wi_e))) );
%scrrecombQE_scr = (1-ref) .* exp(-a_w.*t_w) .* ( 1-exp(-a_i.*(t_i+Wd)) ) .* exp( -a_e .*(t_e-wi_e) ).*exp(-(t_i+Wd-(-1./a_i.*log(((1)-exp(-a_i.*(t_i+Wd))./(a_i.*(t_i+Wd))))))./1.2E7./tauie).*exp(-((-1./a_e.*log(((1)-exp(-a_e.*(t_i+Wd))./(a_e.*(t_i+Wd)))))*(t_i+Wd)/1.2E7/tauip));
QE_scr = (1-ref) .*exp(-2.15E-4.*pass.*a_GaPAbsStaffan).* exp(-a_w.*t_w) .* ( 1-exp(-a_i.*(t_i+Wd)) ) .* exp( -a_e .*(t_e-wi_e) );
QE_b = (1-ref) .* exp(-a_w.*t_w) .* exp(-a_i.*t_i).* ( E./(E.^2-1) .* exp(-a_b.*(t_e+wi_b)) ) .* ( a_b.*Lb - ( ( J.*(cosh((t_b-wi_b)./Lb)-exp(-a_b.*(t_b-wi_b)) ) + sinh((t_b-wi_b)./Lb) + a_b.*Lb.*exp(-a_b.*(t_b-wi_b)) ) ./ (  J.*sinh((t_b-wi_b)./Lb) + cosh((t_b-wi_b)./Lb) ) ) ); 
NQEe=exp(-a_w*t_w).*(1-ref).*a_e.*Le./(a_e.^2*Le^2-1).*( ( Se*Le/De+a_e*Le-exp(-a_e*(t_e-wi_e)).*(Se*Le./De*cosh((t_e-wi_e)/Le)+sinh((t_e-wi_e)/Le)))./(Se*Le/De*sinh((t_e-wi_e)/Le)+cosh((t_e-wi_e)/Le))-a_e.*Le.*exp(-a_e*(t_e-wi_e)))+QEw./(Se.*taue./Le.*sinh(t_e./Le)+cosh(t_e./Le)) ;
NQEscr=(1-ref).*exp(-a_w.*t_w).*(1-exp(-a_i.*(Wd+t_i))).*exp(-a_b.*(t_e-wi_e)) .* exp(-t_i/Li)*exp(-0.5*(Wd+t_i)/1.2E7/taui); %contains attenuation from emitter layer
NQEb=(1-ref).*exp(-2.15E-4.*pass.*a_GaPAbsStaffan).*exp(-a_w.*t_w).*exp(-a_i.*t_i).*a_b*Lb./(a_b.^2*Lb^2-1).*exp(-a_e*t_e-a_b*wi_b).*(a_b*Lb-(Sb*Lb/Db*(cosh((t_b-wi_b)/Lb)-exp(-a_b*(t_b-wi_b)))+sinh((t_b-wi_b)/Lb)+a_b.*Lb.*exp(-a_b*(t_b-wi_b)))./(Sb*Lb/Db*sinh((t_b-wi_b)/Lb)+cosh((t_b-wi_b)/Lb)));

Rback=1;

NQEbR=(1-ref).*exp(-a_w.*t_w).*exp(-a_i.*t_i).*exp(-a_b.*t_b).*a_b*Lb./(a_b.^2*Lb^2-1).*exp(-a_e*t_e-a_b*wi_b).*(a_b*Lb-(Sb*Lb/Db*(cosh((t_b-wi_b)/Lb)-exp(-a_b*(t_b-wi_b)))+sinh((t_b-wi_b)/Lb)+a_b.*Lb.*exp(-a_b*(t_b-wi_b)))./(Sb*Lb/Db*sinh((t_b-wi_b)/Lb)+cosh((t_b-wi_b)/Lb)));
QE_scrR=(1-ref) .* exp(-a_w.*t_w).* exp(-a_b.*2.*t_b).* ( 1-exp(-a_i.*(t_i+Wd)) ) .* exp( -a_e .*(t_e-wi_e) );
QE_eR=(1-ref) .* exp(-a_w.*t_w).*exp(-a_e.*t_e).*exp(-a_i.*2.*t_i).*exp(-a_b.*2.*t_b).* ( S./(S.^2 - 1)) .*  (  ((V+a_e.*Le) - exp(-a_e.*(t_e-wi_e)).*(V.*cosh((t_e-wi_e)./Le)+sinh((t_e-wi_e)./Le) ) ) / ( V.*sinh((t_e-wi_e)./Le)+cosh((t_e-wi_e)./Le) ) - (a_e.*Le.*exp(-a_e.*(t_e-wi_e))) );
QE_wR=(1-ref).* exp(-a_w.*t_w).*exp(-a_e.*2.*t_e).*exp(-a_i.*2.*t_i).*exp(-a_b.*2.*t_b).*a_w.*Lw./(a_w.^2.*Lw.^2-1).*( (a_w.*Lw+Sw.*tauw./Lw.*(1-exp(-a_w.*t_w).*cosh(t_w./Lw))-exp(-a_w.*t_w).*sinh(t_w./Lw))./(Sw.*tauw./Lw.*sinh(t_w./Lw)+cosh(t_w./Lw))-a_w.*Lw.*exp(-a_w.*t_w) );

QE_wR2=(1-ref).* exp(-a_w.*2.*t_w).*exp(-a_e.*2.*t_e).*exp(-a_i.*2.*t_i).*exp(-a_b.*2.*t_b).*a_w.*Lw./(a_w.^2.*Lw.^2-1).*( (a_w.*Lw+Sw.*tauw./Lw.*(1-exp(-a_w.*t_w).*cosh(t_w./Lw))-exp(-a_w.*t_w).*sinh(t_w./Lw))./(Sw.*tauw./Lw.*sinh(t_w./Lw)+cosh(t_w./Lw))-a_w.*Lw.*exp(-a_w.*t_w) );
QE_eR2=(1-ref) .* exp(-a_w.*3.*t_w).*exp(-a_e.*2.*t_e).*exp(-a_i.*2.*t_i).*exp(-a_b.*2.*t_b).* ( S./(S.^2 - 1)) .*  (  ((V+a_e.*Le) - exp(-a_e.*(t_e-wi_e)).*(V.*cosh((t_e-wi_e)./Le)+sinh((t_e-wi_e)./Le) ) ) / ( V.*sinh((t_e-wi_e)./Le)+cosh((t_e-wi_e)./Le) ) - (a_e.*Le.*exp(-a_e.*(t_e-wi_e))) );
QE_scrR2=(1-ref) .* exp(-a_w.*3.*t_w).* exp(-a_b.*2.*t_b).* ( 1-exp(-a_i.*(t_i+Wd)) ) .* exp( -a_e .*3.*(t_e-wi_e) );
NQEbR2=(1-ref).*exp(-a_w.*3.*t_w).*exp(-a_i.*3.*t_i).*exp(-a_b.*2.*t_b).*a_b*Lb./(a_b.^2*Lb^2-1).*exp(-a_e.*3.*t_e-a_b*wi_b).*(a_b*Lb-(Sb*Lb/Db*(cosh((t_b-wi_b)/Lb)-exp(-a_b*(t_b-wi_b)))+sinh((t_b-wi_b)/Lb)+a_b.*Lb.*exp(-a_b*(t_b-wi_b)))./(Sb*Lb/Db*sinh((t_b-wi_b)/Lb)+cosh((t_b-wi_b)/Lb)));


% %before changing absorption coefficients
% NQEe=exp(-a_w*t_w).*(1-ref).*Ap.*Le./(Ap.^2*Le^2-1).*( ( Se*Le/De+Ap*Le-exp(-Ap*(t_e-wi_e)).*(Se*Le./De*cosh((t_e-wi_e)/Le)+sinh((t_e-wi_e)/Le)))./(Se*Le/De*sinh((t_e-wi_e)/Le)+cosh((t_e-wi_e)/Le))-Ap.*Le.*exp(-Ap*(t_e-wi_e)))+QEw./(Se.*taue./Le.*sinh(t_e./Le)+cosh(t_e./Le)) ;
% NQEscr=(1-ref).*exp(-a_w.*t_w).*(1-exp(-a_i.*(Wd+t_i))).*exp(-An.*(t_e-wi_e)); %contains attenuation from emitter layer
% NQEb=(1-ref).*exp(-a_w.*t_w).*exp(-a_i.*t_i).*An*Lb./(An.^2*Lb^2-1).*exp(-Ap*t_e-An*wi_b).*(An*Lb-(Sb*Lb/Db*(cosh((t_b-wi_b)/Lb)-exp(-An*(t_b-wi_b)))+sinh((t_b-wi_b)/Lb)+An.*Lb.*exp(-An*(t_b-wi_b)))./(Sb*Lb/Db*sinh((t_b-wi_b)/Lb)+cosh((t_b-wi_b)/Lb)));

%totals
QE_t =QE_e+QE_w+QE_scr+NQEb;
QE_tR=NQEbR+QE_scrR+QE_eR+QE_wR+QE_t;
QE_tR2=NQEbR2+QE_scrR2+QE_eR2+QE_wR2+QE_tR;


%% Dark Reverse Saturation Currents and 1-Sun JV


Jp=q*De*ni_e^2/(Ne_a*Le)*((sinh((t_e-wi_e)/Le)+Se*Le/De*cosh((t_e-wi_e)/Le))/(cosh((t_e-wi_e)/Le)+Se*Le/De*sinh((t_e-wi_e)/Le)));%*(exp(V./.026)-1)
Jn=q*Db*ni_b^2/(Nb_d*Lb)*((sinh((t_b-wi_b)/Lb)+Sb*Lb/Db*cosh((t_b-wi_b)/Lb))/(cosh((t_b-wi_b)/Lb)+Sb*Lb/Db*sinh((t_b-wi_b)/Lb)));%*(exp(V./.026)-1)
J0=Jp+Jn;

Jsc1sun=sum(F.*dLambda.*QE_t.*Lambda./1240)*1000
V=linspace(0,2,201);
J=J0*(exp(V./.026)-1);
Jlight=Jsc1sun-J*1000;
Voc1sun=.026*log(Jsc1sun/((J0)*1000))
voc=Voc1sun/(Vt);
Rch=Voc1sun/(Jsc1sun/1000);
Rs=0.2;
rs=Rs/Rch;
pmax0=max(Jlight.*V);
FFo=pmax0/(Jsc1sun*Voc1sun);
FFs1sun=FFo*(1-1.1*rs)+rs^2/5.4;%(voc-log(voc+.720))/(voc+1)%
pmaxs=FFs1sun*(Jsc1sun*Voc1sun);
Efficiency1sun=pmaxs/(90)*100;

%% Concentration
[X, Jsc, Voc, FFs, Efficiency]=concFn(F, dLambda, QE_t, Lambda, Jp + Jn);%Concentration Measurements

figure(1)
semilogx(X,FFs)
xlabel('\fontsize{18}Concentration Factor (X)')
ylabel('\fontsize{18}Fill Factor')

figure(2)
semilogx(X,Efficiency)
xlabel('\fontsize{18}Concentration Factor (X)')
ylabel('\fontsize{18}Efficiency (%)')

figure(3)
semilogx(X,Voc)
xlabel('\fontsize{18}Concentration Factor (X)')
ylabel('\fontsize{18}Open Circuit Voltage (V)')

figure(4)
semilogx(X,Jsc)
xlabel('\fontsize{18}Concentration Factor (X)')
ylabel('\fontsize{18}Short Circuit Current Density (mA/cm^2)')

%% Output and Plots
output1=[Lambda QE_w QE_e QE_scr QE_b QE_t]; %Nelson SR outputs

    figure(5)
    plot(Lambda, [QE_w QE_e QE_scr NQEb QE_t])
%    
    title('\fontsize{18}QE Fit')
     xlabel('\fontsize{18}Wavelength (nm)');
    set(gca,'XLim',[300 800],'Layer','top')
    ylabel('\fontsize{18}QE');
    ylim([0 1])
    legend('\fontsize{12}Window','\fontsize{12}Emitter','\fontsize{12}SCR','\fontsize{12}Base', '\fontsize{12}TOTAL') %
    legend('Location',['NorthWest']) %best
   
     figure(6)
    plot(V, [Jlight])
    
    title('\fontsize{18}AM0 JV')
    xlabel('\fontsize{18}Voltage (V)');
    set(gca,'XLim',[0 2],'Layer','top')
    ylabel('\fontsize{18}Current Density (mA/cm^2)');
    ylim([0 5])
 %   legend('\fontsize{12}Model','\fontsize{12}Experimental')
 %   legend('Location',['NorthWest']) %best
   
% end
     