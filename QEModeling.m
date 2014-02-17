% function [ output ] = QEpol(tsrhe, tsrhb,SR09A161, SR10A018C3,Lambda,Fsun, Ref,Aryanadj2,Aw,Ne_a,Naw,Nb_d,wi_b,wi_e,Wd,taue,taub,tauw,mue,mub,muw,Le,Lb,Lw, De,Db,Dw,Se,Sb,Sw,Vt,q,T,C,t_e,tauw,t_b,ni,nib,Vbi,t_i,area,pol)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
clear all
clc
close all
pause(1)

load('~/Dropbox/VM Backups/VM 2013-4-10/Kerestes QE/QEworkspace2.mat'); %absorption, spectra, and reflectance data
%WavelengthEQEIQERSRdata = xlsread('C:\Users\Zac Bittner\Desktop\Kerestes QE\GaAs60XR7.xls'); %file location of spectral response/EQE data
%Wave = WavelengthEQEIQERSRdata(:,1); EQE=WavelengthEQEIQERSRdata(:,3); SR=WavelengthEQEIQERSRdata(:,2);
% load('C:\Users\cxksps\Documents\MATLAB\mats\GaAs_nka3.mat');
% Lnkadata = xlsread('C:\Users\cxksps\Documents\MATLAB\mats\InGaP Epi GenOsc Fit.xls');

%concentration factor
X=1;

%incident solar spectrums W/cm^2/nm for am0, am1.5g, amd1.5d 280-1100nm
spect=1; %spectrum to use [1,2,3] = [AM0,AM1.5G,AM1.5D]
F=Fsun(2:end,spect);
Psun=Psun(spect); %incident solar power W/cm^2 for am0, am1.5g, amd1.5d

%Reflection
ref=Ref;
%ref=0;

T=300 ; %[Kelvin]
Vt=kb*T/q;

%dielectric constants erel_GaAs=13.1, erel_InGaP=11.8
erel_e=erel_InGaP;
erel_i=erel_GaAs; 
erel_b=erel_GaAs; 
% K=13.1; %gaas dielectric constant

% %Polarization: p on n 'pin'  n on p 'nip'
% pol='nip';

%Layer Thicknesses
t_w = 5e-6; %cm
t_e = 5e-5 %cm
t_i = 1.33E-5; %cm 8.73
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
Nb_d = 1e17; %1e17 changed Nd to Nb_d
Nsub=1e18;

%doping (pwindow/pemitter/nbase) [cm-3]
Naw=3.5e18;
Na = 1.2e18;
Ni=1e8;
Nd = 1e17; %1e17
Nsub=1e18;
%mobilities [cm2/Vs]   from Marti paper
%muw=1000; %InGaP window  (minority electrons in ptype)
%mue=250; %pGaAs emitter (minority e- in p type)
%mub=1000; %nGaAs base (minority holes in n type)
%musub=mub
%equations from Marti/Algora1997 
%Electrons.......  as min carriers c=(1+cp/3-cp)=1/3 with no compensation , M=maj carrier concentration
c=1/3;
mumaxn=-6933*c^3+3581*c^2+76*c+9157;
muminn=-2201*c^3+5362*c^2-5092*c+1942;
Nrefn=-4.36e16*c^3+1.43e17*c^2-1.80e17*c+8.49e16;
chi=.506*c^3-.262*c^2+.132*c+.365;
Tn=T/300;
% 
% %Aux equation
% mun=mumin*Tn^-.57+((mumax-mumin)*Tn^-2.3)/(1+(M/(Nref*Tn^-2.4))^(chi*Tn-.146))
%Holes........
mumaxp=400; %[cm2/Vs]
Nrefp=6.25e17; %[cm-3]
%  %Hole aux equation
%mup=mumaxp*Tn^-2.3/(1+M/(Nrefp*Tn^3.8))
Maw=Naw; %maj carrier conc in window
Ma=Na; %maj carrier conc in emitter
Mi=Ni;
Md=Nd; %maj carrier conc in base
Msub=Nsub;

muw=muminn*Tn^(-.57)+((mumaxn-muminn)*Tn^(-2.3))/(1+(Maw/(Nrefn*Tn^(-2.4)))^(chi*Tn^(-.146))); %InGaP window  (minority electrons in ptype)
mue=muminn*Tn.^(-.57)+((mumaxn-muminn)*Tn^(-2.3))./(1+(Ma./(Nrefn*Tn^(-2.4))).^(chi*Tn^(-.146))) %pGaAs emitter (minority e- in p type)
mui=muminn*Tn.^(-.57)+((mumaxn-muminn)*Tn^(-2.3))./(1+(Mi./(Nrefn*Tn^(-2.4))).^(chi*Tn^(-.146))); %pGaAs i-region
mub=mumaxp*(Tn.^(-2.3))./(1+Md./(Nrefp*Tn^3.8)) %nGaAs base (minority holes in n type)
musub=mumaxp*(Tn^(-2.3))/(1+Msub/(Nrefp*Tn^3.8));

%Diffusivity
Dw=Vt*muw;
De=Vt*mue;
Di=Vt.*mui;
Db=Vt*mub;

%Dw = 1; %cm2/s est of AlGaP
%De = 1.086; %cm2/s PC1D InGaP 2e18
%Db = 79.81; %cm2/s PC1D GaAs 1e17
%Dsub = 42.93; %cm2/s PC1D GaAs 1e18

Sw=1E9; %1e8; %front surface InGaP/air interface  %%%% VARIABLE!!!!!!!!!!!!!!!
Se=100; %emitter/window interface recombination  %%%% VARIABLE!!!!!!!!!!!!!!!
Sb=100; %base back surface recombination  %%%% VARIABLE!!!!!!!!!!!!!!!
dataS=[Sw,Se,Sb];

%minority carrier lifetimes
tauw = 2e-13;%s%2E-13
taue = 2.5E-11;%1.4E-11%3.2E-11
taub = 3.5e-9; %s3.5E-11%3.5E-11
tauip = 1.1E-15%1.1E-15 %1.3E-15%1.2E-15;
tauie = 2.6E-11 %2.0E-11%2.1E-11;
taui=tauip;
% tausub1

% %minority carrier lifetimes in the dep region
% tp0=4e-12;       %4, 1e-5%%% VARIABLE!!!!!!!!!!!!!!! ???
% tn0=4e-12;       %4, 1e-5%%% VARIABLE!!!!!!!!!!!!!!! ???
% C=1/sqrt(tn0*tp0);
% C=1e9; %from Algora reference 1e5, Hovel C=1e9 s^-1

Lw=sqrt(Dw.*tauw)%cm
Le=sqrt(De.*taue)%cm
Li=2e-4; %cm 
Lb=sqrt(Db.*taub)%cm
dataL=[Lw,Le,Lb];

%Temp and Concentration dependent ni for GaAs
Eg_GaAs=1.519-5.405*10^-4*T^2/(204+T); %in eV, for intrinsic GaAs  
Nc_GaAs=8.63e13*T^(3/2)*(1-1.93e-4*T-4.19e-8*T^2);
Nv_GaAs=1.83e15*T^(3/2);
dEg_GaAs=2e-11*Ne_a.^(.5); %in eV, p-type gap narrowing from Marti
dEgn_GaAs=.125/erel_e*abs(1-Nb_d./Ncr).^(1/3); %eV (intrinsic (n-type?) gap narrowing from Marti)
ni_GaAs=((Nc_GaAs*Nv_GaAs)^.5)*exp((-Eg_GaAs-dEgn_GaAs)./(2*Vt)); %ni in base

ni_w = 1500;
ni_e = ni_GaAs; %Ga(x)In(1-x)P, x=51.5% Eg=1.85
ni_b = ni_GaAs;

%built in junction voltage
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


%Absorption Coefficients
a_w=a_InGaPEPI; %absorption in windowchanged Aw to a_w
a_e=a_iGaAsEPI; %absorption in emitterchanged An to a_e
a_i=a_iGaAsEPI; %absorption in i-regionGaAs changed Ai to a_i
a_b=a_iGaAsEPI; %absorption in baseGaAs changed Ap to a_b

%subscript w=window (a in Hovel) e=emitter (g in Hovel), b=base (p Hovel)
%absorption coefficients taken from Aguinaldo data p=pGaAs, n=nGaAs, and
%measured GaAs data from Michael Slocum

% if pol=='nip'
% %     Bhold=tsrhe;tsrhe=tsrhb;tsrhb=Bhold;
%     Bhold=Ne_a;Ne_a=Nb_d;Nb_d=Bhold;
% %     Bhold=Naw;Naw=Ndw;Ndw=Bhold; window layer
%     Bhold=wi_b; wi_b=wi_e; wi_b=wi_e;
% %     Bhold=taue; taue=taub; taub=Bhold;
% %     Bhold=mue; mue=mub; mub=Bhold;
%     Bhold=Le;Le=Lb; Lb=Bhold;      
%     Bhold=De; De=Db; Db=Bhold;
%     Bhold=Se; Se=Sb; Sb=Bhold;
%     Bhold=a_e; a_e=a_b; a_b=Bhold;
% %     Bhold=t_e; t_e=tw; tw=Bhold;
%     %nib : need equation for band gap narrowing of p-type GaAs nib is Eg
%     %narrowing effects for n-type GaAs
% elseif pol=='pin'
%     ZZZ=1; %no changes made to variables in MaterialParams.m
% else
%     disp(sprintf('Error: pin or nip polarity'))
% end

% Constants for QE fitting in window layer
Q=a_w*Lw; %changed B to a_w
Y=-t_w/Lw; %changed D to t_w
Z=-Sw*tauw/Lw; %chaged tw to tauw

% Constants for QE fitting in emitter layer
S=a_e.*Le;
T=-t_e/Le; %changed d to t_e
U=-Se*taue/Le; %changed te to taue
V = Se*Le/De;

% Constants for QE fitting in base layer
E=a_b.*Lb;
Fa=t_b/Lb; %changed Tb to t_b
G=Sb.*taub./Lb; %changed tb to taub
J = Sb.*Lb./Db;




QEwAlgora=-(1-ref).*exp(-a_w.*t_w).*(Q./(Q.^2-1).*(Q-1+( (Z-Q).*exp(-Q.*Y)+(1-Z).*exp(-Y))./(cosh(Y)+Z*sinh(Y))));
QEeAlgora=-(1-ref).*exp(-a_w.*t_w-a_e.*t_e).*(S./(S.^2-1).*(S-1+( (U-S).*exp(-S.*T)+(1-U).*exp(-T))./(cosh(T)+U*sinh(T))));
QEbAlgora=(1-ref).*exp(-a_w.*t_w - a_e.*t_e - a_b.*Wd - a_i.*t_i).*(E./(E.^2-1).*(E-1+( (G-E).*exp(-E.*Fa)+(1-G).*exp(-Fa))./(cosh(Fa)+G*sinh(Fa))));

QEw=(1-ref).*a_w.*Lw./(a_w.^2.*Lw.^2-1).*( (a_w.*Lw+Sw.*tauw./Lw.*(1-exp(-a_w.*t_w).*cosh(t_w./Lw))-exp(-a_w.*t_w).*sinh(t_w./Lw))./(Sw.*tauw./Lw.*sinh(t_w./Lw)+cosh(t_w./Lw))-a_w.*Lw.*exp(-a_w.*t_w) );
QEe=(1-ref).*exp(-a_w.*t_w).*a_e.*Le./(a_e.^2.*Le.^2-1).*((a_e.*Le+Se.*taue./Le.*(1-exp(-a_e.*t_e).*cosh(t_e./Le))-exp(-a_e.*t_e).*sinh(t_e./Le))./(Se.*taue./Le.*sinh(t_e./Le)+cosh(t_e./Le))-a_e.*Le.*exp(-a_e.*t_e))+QEw./(Se.*taue./Le.*sinh(t_e./Le)+cosh(t_e./Le))*exp(-(Wd+t_i)/1.2E7/tauie);
QEb=(1-ref).*exp(-a_w.*t_w).*a_b.*Lb./(a_b.*Lb+1).*exp(-a_e.*t_e).*exp(-a_i.*Wd); %from Hovel
QEscr=(1-ref).*exp(-a_w.*t_w).*(1-exp(-a_i.*Wd)).*exp(-a_e.*(t_e)); %from Hovel (approximates t_e>>wi_e) nearly the same as Nelson
% %before changing absorption coefficients
% QEw=(1-ref).*a_w.*Lw./(a_w.^2.*Lw.^2-1).*( (a_w.*Lw+Sw.*tauw./Lw.*(1-exp(-a_w.*t_w).*cosh(t_w./Lw))-exp(-a_w.*t_w).*sinh(t_w./Lw))./(Sw.*tauw./Lw.*sinh(t_w./Lw)+cosh(t_w./Lw))-a_w.*Lw.*exp(-a_w.*t_w) );
% QEe=(1-ref).*exp(-a_w.*t_w).*Ap.*Le./(Ap.^2.*Le.^2-1).*((Ap.*Le+Se.*taue./Le.*(1-exp(-Ap.*t_e).*cosh(t_e./Le))-exp(-Ap.*t_e).*sinh(t_e./Le))./(Se.*taue./Le.*sinh(t_e./Le)+cosh(t_e./Le))-Ap.*Le.*exp(-Ap.*t_e))+QEw./(Se.*taue./Le.*sinh(t_e./Le)+cosh(t_e./Le)) ;
% QEb=(1-ref).*exp(-a_w.*t_w).*An.*Lb./(An.*Lb+1).*exp(-An.*t_e).*exp(-An.*Wd); %from Hovel
% QEscr=(1-ref).*exp(-a_w.*t_w).*(1-exp(-An.*Wd)).*exp(-An.*(t_e)); %from Hovel (approximates t_e>>wi_e) nearly the same as Nelson

%average absorption distance in i-region
%xavg=zeros(941,1);
%xavg=(-1./a_e.*log((1-exp(-a_e.*(t_i+Wd))./(a_e.*(t_i+Wd)))));
t_trans=(Wd+t_i)/1.2E7
%from Nelson:
QE_w=(1-ref).*a_w.*Lw./(a_w.^2.*Lw.^2-1).*( (a_w.*Lw+Sw.*tauw./Lw.*(1-exp(-a_w.*t_w).*cosh(t_w./Lw))-exp(-a_w.*t_w).*sinh(t_w./Lw))./(Sw.*tauw./Lw.*sinh(t_w./Lw)+cosh(t_w./Lw))-a_w.*Lw.*exp(-a_w.*t_w) )*exp(-(Wd+t_i)/1.2E7/tauie);
QE_e = (1-ref) .* exp(-a_w.*t_w) .* ( S./(S.^2 - 1)) .*  (  ((V+a_e.*Le) - exp(-a_e.*(t_e-wi_e)).*(V.*cosh((t_e-wi_e)./Le)+sinh((t_e-wi_e)./Le) ) ) / ( V.*sinh((t_e-wi_e)./Le)+cosh((t_e-wi_e)./Le) ) - (a_e.*Le.*exp(-a_e.*(t_e-wi_e))) )*exp(-(Wd+t_i)/1.2E7/tauie);
%scrrecombQE_scr = (1-ref) .* exp(-a_w.*t_w) .* ( 1-exp(-a_i.*(t_i+Wd)) ) .* exp( -a_e .*(t_e-wi_e) ).*exp(-(t_i+Wd-(-1./a_i.*log(((1)-exp(-a_i.*(t_i+Wd))./(a_i.*(t_i+Wd))))))./1.2E7./tauie).*exp(-((-1./a_e.*log(((1)-exp(-a_e.*(t_i+Wd))./(a_e.*(t_i+Wd)))))*(t_i+Wd)/1.2E7/tauip));
QE_scr = (1-ref) .* exp(-a_w.*t_w) .* ( 1-exp(-a_i.*(t_i+Wd)) ) .* exp( -a_e .*(t_e-wi_e) );
QE_b = (1-ref) .* exp(-a_w.*t_w) .* exp(-a_i.*t_i).* ( E./(E.^2-1) .* exp(-a_b.*(t_e+wi_b)) ) .* ( a_b.*Lb - ( ( J.*(cosh((t_b-wi_b)./Lb)-exp(-a_b.*(t_b-wi_b)) ) + sinh((t_b-wi_b)./Lb) + a_b.*Lb.*exp(-a_b.*(t_b-wi_b)) ) ./ (  J.*sinh((t_b-wi_b)./Lb) + cosh((t_b-wi_b)./Lb) ) ) ); 
NQEe=exp(-a_w*t_w).*(1-ref).*a_e.*Le./(a_e.^2*Le^2-1).*( ( Se*Le/De+a_e*Le-exp(-a_e*(t_e-wi_e)).*(Se*Le./De*cosh((t_e-wi_e)/Le)+sinh((t_e-wi_e)/Le)))./(Se*Le/De*sinh((t_e-wi_e)/Le)+cosh((t_e-wi_e)/Le))-a_e.*Le.*exp(-a_e*(t_e-wi_e)))+QEw./(Se.*taue./Le.*sinh(t_e./Le)+cosh(t_e./Le))*exp(-(Wd+t_i)/1.2E7/tauie) ;
NQEscr=(1-ref).*exp(-a_w.*t_w).*(1-exp(-a_i.*(Wd+t_i))).*exp(-a_b.*(t_e-wi_e)) .* exp(-t_i/Li)*exp(-0.5*(Wd+t_i)/1.2E7/taui); %contains attenuation from emitter layer
NQEb=(1-ref).*exp(-a_w.*t_w).*exp(-a_i.*t_i).*a_b*Lb./(a_b.^2*Lb^2-1).*exp(-a_e*t_e-a_b*wi_b).*(a_b*Lb-(Sb*Lb/Db*(cosh((t_b-wi_b)/Lb)-exp(-a_b*(t_b-wi_b)))+sinh((t_b-wi_b)/Lb)+a_b.*Lb.*exp(-a_b*(t_b-wi_b)))./(Sb*Lb/Db*sinh((t_b-wi_b)/Lb)+cosh((t_b-wi_b)/Lb)))*exp(-(Wd+t_i)/1.2E7/tauip);

% %before changing absorption coefficients
% NQEe=exp(-a_w*t_w).*(1-ref).*Ap.*Le./(Ap.^2*Le^2-1).*( ( Se*Le/De+Ap*Le-exp(-Ap*(t_e-wi_e)).*(Se*Le./De*cosh((t_e-wi_e)/Le)+sinh((t_e-wi_e)/Le)))./(Se*Le/De*sinh((t_e-wi_e)/Le)+cosh((t_e-wi_e)/Le))-Ap.*Le.*exp(-Ap*(t_e-wi_e)))+QEw./(Se.*taue./Le.*sinh(t_e./Le)+cosh(t_e./Le)) ;
% NQEscr=(1-ref).*exp(-a_w.*t_w).*(1-exp(-a_i.*(Wd+t_i))).*exp(-An.*(t_e-wi_e)); %contains attenuation from emitter layer
% NQEb=(1-ref).*exp(-a_w.*t_w).*exp(-a_i.*t_i).*An*Lb./(An.^2*Lb^2-1).*exp(-Ap*t_e-An*wi_b).*(An*Lb-(Sb*Lb/Db*(cosh((t_b-wi_b)/Lb)-exp(-An*(t_b-wi_b)))+sinh((t_b-wi_b)/Lb)+An.*Lb.*exp(-An*(t_b-wi_b)))./(Sb*Lb/Db*sinh((t_b-wi_b)/Lb)+cosh((t_b-wi_b)/Lb)));

%totals
QEtAlgora=QEeAlgora+QEwAlgora+QEbAlgora+NQEscr;
QE_t =QE_e+QE_w+QE_scr+NQEb;
NQEt=NQEe+NQEb+NQEscr;
QEt=QEe+QEb+QEscr; %QEw+

SRtAlgora=QEtAlgora.*Lambda./1240;
SR=QEt.*Lambda./1240; %make sure this is correct
NSR=NQEt.*Lambda./1240;

%individual SR contributions
AlgoraSRw=QEwAlgora.*Lambda./1240;
AlgoraSRe=QEeAlgora.*Lambda./1240;
AlgoraSRb=QEbAlgora.*Lambda./1240;
AlgoraSRscr=NQEscr.*Lambda./1240;

SRw=QEw.*Lambda./1240;
SRe=QEe.*Lambda./1240;
SRb=QEb.*Lambda./1240;
SRscr=QEscr.*Lambda./1240;
NSRb=NQEb.*Lambda./1240;
NSRscr=NQEscr.*Lambda./1240;
NSRe=NQEe.*Lambda./1240;

% %Sets Hovel Model equal to Nelson Model?????????
% QEb=NQEb;
% QEscr=NQEscr;
% QEt=NQEt;

%plot(Lambda(641:841),[QE09A161(641:841) NQEt(641:841)]); figure;
%plot(Lambda(641:841), QDT(641:841,9)/100/normalization, 'magenta'); hold on;
%plot(Lambda(641:841), BaselineT(641:841,9)/100/normalization, 'blue'); hold on;
%plot(Lambda(641:841), [NQEb(641:841) NQEscr(641:841) QEe(641:841) NQEt(641:841)])
%plot(Lambda, [NSRscr NSRb SRe NSRscr+NSRb+SRe SR09A161])
%plot(Lambda,[SRe SRb SRscr SR SR09A161]);
%plot(Lambda, [SR09A161 SR])
%Photocurrent summation
%plot(Lambda(200:end), [NSR(200:end) SRw(200:end) NSRb(200:end) SRe(200:end) NSRscr(200:end) SR09A161(200:end)])
%plot(Lambda, [SR10A018C3 NSR]);
%plot(Lambda, [SR10A018C3 NSR]);

%% For calculating light and dark current values
% Lambda2=circshift(Lambda,1);
% Lambda2(1)=[];
% dLambda=Lambda(2:end)-Lambda2;
% jphi=X*sum(F.*dLambda.*(NSR(2:end)));% Integrated Current Density
% jphi09A161=X*sum(F.*dLambda.*(SR09A161(2:end)));% Integrated Current Density
% jphi10A018C3=X*sum(F.*dLambda.*(SR10A018C3(2:end)));% Integrated Current Density
% jphiam0=X*sum(F.*dLambda.*(NSR(2:end)));% Integrated Current Density
% % Integrated Current Density of individual current components
% jphe=X*sum(F.*dLambda.*(NSRe(2:end)));
% jphb=X*sum(F.*dLambda.*(NSRb(2:end)));
% jphscr=X*sum(F.*dLambda.*(NSRscr(2:end)));


 %Ploted output: SRw NSRe NSRb NSRscr NSR(From Nelson)
 output1=[Lambda QE_w QE_e QE_scr NQEb QE_t]; %Nelson SR outputs
%  output1=[Lambda QEwAlgora QEeAlgora QEscr QEbAlgora SRwAlgora SReAlgora SRbAlgora SRscr SRAlgora]; %Algora SR outputs
%  output1=[Lambda QEw QEe QEscr QEb SRw SRe SRb SRscr SR]; %Hovel SR outputs

    figure(1)
    plot(Lambda, [QE_w QE_e QE_scr NQEb QE_t])
    
    title('\fontsize{18}QE Fit')
    xlabel('\fontsize{18}Wavelength (nm)');
    set(gca,'XLim',[300 1000],'Layer','top')
    ylabel('\fontsize{18}QE');
    ylim([0 .8])
    legend('\fontsize{12}Window','\fontsize{12}Emitter','\fontsize{12}SCR','\fontsize{12}Base', '\fontsize{12}TOTAL' ) %
    legend('Location',['NorthWest']) %best
   
% end
     