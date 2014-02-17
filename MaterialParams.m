clear all
clc
close all
pause(1)
%%%%%%%%%%%%%%%%%%%%INPUT PARAMS%%%%%%%%%%%%%%%%%%%%%%%

%e=emitter b=base w=window  [g=emitter (pGaAs)  p=base (nGaAs) a=window
%(pInGaP) Hovel]
%data in variable SR10A018C3 is good SR data, only active area used,
%integrated SR matches jsc value

load('C:\Users\Zac Bittner\Desktop\Kerestes QE\QEworkspace.mat');
load('C:\Users\Zac Bittner\Desktop\Kerestes QE\QEworkspace2.mat');
load('C:\Users\Zac Bittner\Desktop\Kerestes QE\WavelengthRefdata.txt'); %file location of reflection data
% load('C:\Users\cxksps\Documents\MATLAB\HarrisQEModel\Data\WavelengthSREQEdata.txt'); %file location of spectral response/EQE data
WavelengthSREQEdata = xlsread('C:\Users\Zac Bittner\Desktop\Kerestes QE\UCLA2012-3-6.xlsx'); %file location of spectral response/EQE data
SR=WavelengthSREQEdata(:,2);
% nmnka_InAs_Palik = xlsread('C:\Users\cxksps\Documents\MATLAB\mats\nmnka_InAs-Palik.xls'); %file location of spectral response/EQE data

gg=1;
%for L=.1:.1:5
% gg
% L=1;
%what data should QE function return: 'IV' for IV data and power loss, 'SR' for spectral
%responses, 'dark' for IV curve and dark currents 
returnVal='QE';  
uu=returnVal;
%Matrix8=[];
%variable=[];
%for tf=.5e-4:.5e-4:10e-4
   
%constants [si units + cm]
% ref=.34; %reflectivity of surface
% nk = padarray(1, [0 (length(n_GaAs)-1)],'symmetric','post'); %probability that carriers are collected
%CALC REFLECTANCE

spect=1; %spectrum to use [1,2,3] = [AM0,AM1.5G,AM1.5D]
T=300 ; %[Kelvin]
q=1.602e-19;
kb=1.381e-23; %boltzmann
Vt=kb*T/q;
erel=12.4*(1+1.2e-4); %relative dielectric constant GaAs (from Marti)
erel=13.1;
K=13.1; %gaas dielectric constant
e0=8.85e-14; %free space permittivity

Ncr=1.6e24*(.066/(1.4*erel))^3;
%doping (pwindow/pemitter/nbase) [cm-3]
Naw=3.5e18;
Na = 1.2e18;
Ni=1e15;
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
mue=muminn*Tn.^(-.57)+((mumaxn-muminn)*Tn^(-2.3))./(1+(Ma./(Nrefn*Tn^(-2.4))).^(chi*Tn^(-.146))); %pGaAs emitter (minority e- in p type)
mui=muminn*Tn.^(-.57)+((mumaxn-muminn)*Tn^(-2.3))./(1+(Mi./(Nrefn*Tn^(-2.4))).^(chi*Tn^(-.146))); %pGaAs i-region
mub=mumaxp*(Tn.^(-2.3))./(1+Md./(Nrefp*Tn^3.8)); %nGaAs base (minority holes in n type)
musub=mumaxp*(Tn^(-2.3))/(1+Msub/(Nrefp*Tn^3.8));

%mue=mub
%Surface recomb
Sw=1e8; %1e8; %front surface InGaP/air interface  %%%% VARIABLE!!!!!!!!!!!!!!!
Se=2e5; %2e5 emitter/window interface recombination  %%%% VARIABLE!!!!!!!!!!!!!!!
% Sei=; %emitter/i-region interface recombination
% Sib=; %i-region/base interface recombination
Sb=1e6; %base back surface recombination  %%%% VARIABLE!!!!!!!!!!!!!!!
dataS=[Sw,Se,Sb];

%lifetimes, C=1/sqrt(n0p0) [s]
%minority e- (p-window)
tradw=1/(2e-10*Naw);
taugw=1/(1.10e-30*Naw^2);
tsrhw=1.3e-13; %1e-13 5e-8;             %%%% VARIABLE!!!!!!!!!!!!!!!
%p-window sum
tw=1/(1/tradw+1/taugw+1/tsrhw);%1.5e-13

%Marti trad^^^
%minority e- (emitter)
trade=1./(1.7e-10.*Na);
tauge=1./(1.10e-30*Na.^2);
tsrhe=2.3e-11; %5e-11(60X) 5e-8(10X);             %%%% VARIABLE!!!!!!!!!!!!!!!
%p-emitter sum
te=1./(1./trade+1./tauge+1./tsrhe);

%minority carrier lifetimes in the dep region
tp0=4e-12;       %4, 1e-5%%% VARIABLE!!!!!!!!!!!!!!! ???
tn0=4e-12;       %4, 1e-5%%% VARIABLE!!!!!!!!!!!!!!! ???
C=1/sqrt(tn0*tp0);
C=1e9; %from Algora reference 1e5, Hovel C=1e9 s^-1

%based on dopings (Marti)
%minority holes (base) 
tradb=1./(1.7e-10*Nd);
taugb=1./(4.72e-30*Nd.^2);
tsrhb=1e-3; %4e-8;             %%%% VARIABLE!!!!!!!!!!!!!!!
%tradbm=10^(.00455*log10(Nd)^3+.608*log10(Nd)^2-26.706*log10(Nd)+249.89)
%tradb=tradbm
%n-base sum
tb=1./(1./tradb+1./taugb+1./tsrhb);

% te=te/10;
%fitting 
%tb=te
%tb=1e-9
%tw=tw;
%tb=3*tb;
% tb=10e-8
% te=10e-8
% tw=10e-8
%tb=2e-8;
%mub=350;
%te=1.5e-9;
%tw=1e-11; window lifetime has very little effect on SR
%muw=1
%mue=150;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%   STRUCTURE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%window thickness [cm]
D=50e-7;       %%%% VARIABLE!!!!!!!!!!!!!!!
%emitter thickness [cm]
d=500e-7;       %%%% VARIABLE!!!!!!!!!!!!!!!
%i-region thickness [cm]
%0X: 100e-7; 5X: 133.25e-7; 10X: 205e-7; 20X: 335e-7; 40X: 604e-7;
%60X:873e-7; 100X: 1411e-7;
wi=873e-7; % 100nm -width of intrinsic layer, approximately adds to dep layer        %%%% VARIABLE!!!!!!!!!!!!!!!
%base thickness [cm]
Tb=2e-4;       %%%% VARIABLE!!!!!!!!!!!!!!!
%GaAs thickness [cm] (H-D) in Hovel
%G=Tb+d;
%substrate thickness [cm]
Tsub=350e-4;       %%%% VARIABLE!!!!!!!!!!!!!!!
rw=2*2.54; %radius of 2inch wafer [cm]
areaw=pi*rw^2; %area of wafer for substrate resistivity [cm2]

%calculated
%tb=5e-8
%mub=600
%diffusion Co 
Dw=Vt*muw;
De=Vt*mue;
Di=Vt.*mui;
Db=Vt*mub;

%diffusion lengths
Lw=sqrt(Dw.*tw);
Le=sqrt(De.*te);
Li=sqrt(Di.*tn0);
Lb=sqrt(Db.*tb);

%estimates
%Lw=1e-4;
%Le=5e-4;
%Lb=4e-4;
%BandGap dependence on T
Eg=1.519-5.405*10^-4*T^2/(204+T); %in eV, for intrinsic GaAs  
Nc=8.63e13*T^(3/2)*(1-1.93e-4*T-4.19e-8*T^2);
Nv=1.83e15*T^(3/2);
dEg=2e-11*Na.^(.5); %in eV, gap narrowing from Marti
dEgn=.125/erel*abs(1-Nd./Ncr).^(1/3); %eV (intrinsic gap narrowing from Marti)
%What is equation for Eg narrowing of p-type GaAs???????????????????
dEgn=0;
ni=((Nc*Nv)^.5).*exp((-Eg-dEg)./(2*Vt));
nib=((Nc*Nv)^.5)*exp((-Eg-dEgn)./(2*Vt)); %ni in base

%intrinsic concentration (with T and band gap narrowing)
%ni0=2.1e6;
%Mp
%Mn
%ni=sqrt(ni0^2*exp(-dEg/(kb*T)));

%built in junction voltage
Vbi=Vt*log(Nd.*Na./ni.^2);
Vbin=Vt*log(Nd/nib);
Vbip=Vt*log(Na/ni);
Vbisep=Vt*( log(Nd./nib)+log(Na./ni) );
Vbi=Vbisep;
%depletion width 0V
Wdi=(2*K*e0/q*(Na+Nd)./(Na.*Nd).*(Vbi-.95)).^.5;
WdiBart=(wi.^2 + (2.*K.*e0./q) .* (Na+Nd)./(Na.*Nd) .* (Vbi-0) ).^.5 - wi; %from Bart
% Wdi=WdiBart;
Vpt = wi.^2 * 1e15* q ./ (2*K*e0);

wp=1./Na.*sqrt(2*K*e0*Vbi./(q*(1./Na+1./Nd)));
wn=1./Nd.*sqrt(2*K*e0*Vbi./(q*(1./Na+1./Nd)));
Wd=wp+wn;
%Wd=Wd/100
%wi=1e-7
%dark current stuff
sL=3;       %%%% VARIABLE!!!!!!!!!!!!!!! WHAT DOES THIS DO???????

% grid stuff       %%%% VARIABLES!!!!!!!!!!!!!!!
L=1;  %length parallel to busbar, perpendicular to fingers, "length" is also function name
width=1;  %width along finger direction
Lper=2*(L+width);
Lper=.45
area=L*width;
rhosilver=1.59e-6; %silver resistivity [ohm-cm]
rhogold=2.27e-6; %gold resistivity [ohm-cm], recently measured ~2.3 ohm-cm
tf=2e-4; %finger height [cm]       %%%% VARIABLE!!!!!!!!!!!!!!!
wf=8e-4; %finger width [cm]       %%%% VARIABLE!!!!!!!!!!!!!!!
rhom=rhogold; %metal resistivity [ohm-cm]
rhoe=1./(q*mue.*(Na)); %emitter resistivity [ohm-cm], decreases with higher J
rhob=1./(q*mub.*Nd); %base resistivity [ohm-cm]
rhoi=3.3e8; %intrinsic resistivity [ohm-cm]       %%%% VARIABLE!!!!!!!!!!!!!!!
rhosub=1/(q*musub*Nsub); %substrate resistivity [ohm-cm]
rhoc=5e-5; %specific contact resistance [ohm-cm^2]        %%%% VARIABLE!!!!!!!!!!!!!!!
Rshe=rhoe./d;%emitter sheet res=rhoe/(emitter thickness d) [ohm/sq]
Rshe=400;
Rshm=rhom/tf; %metal sheet res= rhom/(finger thickness t) [ohm/sq]
%area;
Rshe=200;


%Matrix4=[te tb mue mub De Db Le Lb Na Nd];
%Matrix4=transpose(Matrix4);
%Matrix6=[jphb jphscr jph]*1000;

inputs=[tsrhe;	tsrhb;	areaw;	Tsub;	rhosub;	rhob;	sL;	Lper;	Na;	Naw;	Nd;	wn;	wp;	Wd;	te;	tb;	tw;	mue;	mub;	muw;	Le;	Lb;	Lw;	De;	Db;	Dw;	Se;	Sb;	Sw;	Vt;	q;	T;	C;	d;	D;	Tb;	ni;	nib;	Vbi;	wi;	area;	L;	width;	rhoc;	Rshe;	Rshm;	tf;	wf;	area;	L;	width;	rhoc;	Rshe;	Rshm;	tf;	wf];

variable=QEpol(uu , tsrhe , tsrhb, areaw,Tsub,rhosub,rhob,sL,Lper,SR,Lambda,Fsun, V,Ref,Aryanadj2,B,Na,Naw,Nd,wn,wp,Wd,te,tb,tw,mue,mub,muw,Le,Lb,Lw, De,Db,Dw,Se,Sb,Sw,Vt,q,T,C,d,D,Tb,ni,nib,Vbi,wi,area,L,width,rhoc,Rshe,Rshm,tf,wf,'pin');

% variable=from QE.m:
% SR/QE output1=[Lambda QEw NQEe NQEscr NQEb (6)SRw (7)NSRe (8)NSRb (9)NSRscr (10)NSR];
% QE/SR with 'N' denotes equation from Nelson without 'H' denotes Hovel model
% IV output2=[Effi FFi Eff FF Pmax jmax vmax voc jsc area J01e J01b J02dep(1) J02perim sp wb ps psb pc pbc pe prf prb prbase prsub numf]; %IV performance params and grid
% dark output3=[V Jt Jtpos J1e J1b J2dep J2perim power ]; %IV curve, dark current

if returnVal=='QE'
    figure(1)
    plot(Lambda, [variable(:,6).*1240./Lambda variable(:,7).*1240./Lambda variable(:,8).*1240./Lambda variable(:,9).*1240./Lambda variable(:,10).*1240./Lambda], WavelengthSREQEdata(:,1),WavelengthSREQEdata(:,2),'*k')
    title('\fontsize{18}QE Fit')
    xlabel('\fontsize{18}Wavelength (nm)');
    set(gca,'XLim',[300 1000],'Layer','top')
    ylabel('\fontsize{18}QE');
    YLim([0 1])
    legend('\fontsize{12}Window (Hovel)', '\fontsize{12}Emitter (Nelson)', '\fontsize{12}Base (Nelson)', '\fontsize{12}SCR (Nelson)', '\fontsize{12}TOTAL (Nelson)', '\fontsize{12}Experimental Data' ) %
    legend('Location',['NorthWest']) %best
    
%     figure(2)
%     % plot(Lambda, [SRw SRe SRb SRscr SR], Lambda, [SRw NSRe NSRb NSRscr NSR], WavelengthSREQEdata((3:83),1),WavelengthSREQEdata((3:83),2),'*k');
%     plot(Lambda, [variable(:,6) variable(:,7) variable(:,8) variable(:,9) variable(:,10)], WavelengthSREQEdata(:,1),WavelengthSREQEdata(:,2),'*k')
%     title('\fontsize{18}SR Fit')
%     xlabel('\fontsize{18}Wavelength (nm)');
%     set(gca,'XLim',[300 1000],'Layer','top')
%     ylabel('\fontsize{18}Spectral Response (A/W)');
%     YLim([0 1])
%     legend('\fontsize{12}Window (Hovel)', '\fontsize{12}Emitter (Nelson)', '\fontsize{12}Base (Nelson)', '\fontsize{12}SCR (Nelson)', '\fontsize{12}TOTAL (Nelson)', '\fontsize{12}Experimental Data' ) %
%     legend('Location',['best'])
    
    %[Window,Emitter,Base,SCR];SRV(cm/s);Lifetime(s);Diffusion Length (cm), Thickeness(cm);SRdata(A/W)
    data = [Sw,Se,Sb,0,0,0; tw,te,tb,tp0,0,0; Lw,Le,Lb,0,0,0; D,d,Tb,Wd,0,0; Tsub,0,0,0,0,0; variable(:,6),variable(:,7),variable(:,8),variable(:,9),Lambda, variable(:,10)];
    
elseif returnVal=='IV'

    
    d1={'area(cm)';'voc(V)';'jsc(A/cm2)';'FF(ratio)';'Eff(ratio)';'vmax(V)';'jmax(A/cm2)';'Pmax(W)'; 'Effi';'FFi';'J01e(A/cm2)';'J01b(A/cm2)';'J02dep(1)(A/cm2)';'J02perim(A/cm2)'};
    data1=[area;voc;jsc;FF;Eff;vmax;jmax;Pmax; Effi;FFi;J01e;J01b;J02dep(1);J02perim];
    dataB =[0;0;0;0;0;0;0;0;0;0;0;0;0;0];
    d2={'sp';'wb';'ps';'psb'; 'pc';'pbc';'pe';'prf'; 'prb';'prbase';'prsub';'numf';'0';'0'};
    data2=[sp;wb;ps;psb; pc;pbc;pe;prf; prb;prbase;prsub;numf;0;0]; %IV performance params and grid
    data = [data1,data2];
    datad = [d1,d2];
    
    variable=QEpol('dark' , tsrhe , tsrhb, areaw,Tsub,rhosub,rhob,sL,Lper,SR09A161, SR10A018C3,Lambda,Fsun, V,Ref,Aryanadj2,B,Na,Naw,Nd,wn,wp,Wd,te,tb,tw,mue,mub,muw,Le,Lb,Lw, De,Db,Dw,Se,Sb,Sw,Vt,q,T,C,d,D,Tb,ni,nib,Vbi,wi,area,L,width,rhoc,Rshe,Rshm,tf,wf,'pin');
    data3=[V, Jt, Jtpos, J1e, J1b, J2dep, J2perim, power ];
    d3={'V', 'Jt', 'Jtpos', 'J1e', 'J1b', 'J2dep', 'J2perim', 'power'};
    
        figure(1)
        plot(variable(:,1), variable(:,2))
        title('\fontsize{18}Simulated LIV')
        xlabel('\fontsize{18}Voltage (V)');
        set(gca,'XLim',[0 voc],'Layer','top')
        ylabel('\fontsize{18}Current (A)');
        YLim([0 jsc])
        legend('\fontsize{12}LIV)') %
        legend('Location','best')
        
        figure(2)
        semilogy(variable(:,1), [variable(:,4),variable(:,5),variable(:,6),variable(:,7)])
        legend('\fontsize{12}J1e','\fontsize{12}J1b','\fontsize{12}J2dep','\fontsize{12}J2perim')
pause(1)
        fprintf('Voc:%1.2f V\nIsc:%1.2f mA/cm2\nFF:%1.1f%%\nEff:%1.1f%%',voc, jsc*1000,FF*100,Eff*100)

end

gg=gg+1;
%end
%plot(1e17:1e17:1e18,Matrix8(:,1),'green' );

%end
