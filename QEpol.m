function [ output ] = QEpol(uu, tsrhe, tsrhb,areaw,Tsub,rhosub,rhob,sL,Lper,SR09A161, Lambda,Fsun, V,Ref,Aryanadj2,B,Na,Naw,Nd,wn,wp,Wd,te,tb,tw,mue,mub,muw,Le,Lb,Lw, De,Db,Dw,Se,Sb,Sw,Vt,q,T,C,d,D,Tb,ni,nib,Vbi,wi,area,L,width,rhoc,Rshe,Rshm,tf,wf,pol)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%concentration factor
X=1;
%grid shadowing
ps=0;

%incident solar spectrums W/cm^2/nm for am0, am1.5g, amd1.5d 280-1100nm
F1=Fsun(2:end,1); %AM0
F2=Fsun(2:end,2) ;%AM1.5G
F3=Fsun(2:end,3); %AM1.5D
F=F1;
%incident solar power W/cm^2 for am0, am1.5g, amd1.5d, choose based on
%incident spectrum above
Pam0=.136;
Pam15G=.1;
Pam15D=.09;
Psun=Pam0; 
ref=Ref;
%ref=0;

if pol=='nip'
    Bhold=tsrhe;tsrhe=tsrhb;tsrhb=Bhold;
    Bhold=Na;Na=Nd;Nd=Bhold;
%     Bhold=Naw;Naw=Ndw;Ndw=Bhold; window layer
    Bhold=wn; wn=wp; wn=wp;
    Bhold=te; te=tb; tb=Bhold;
    Bhold=mue; mue=mub; mub=Bhold;
    Bhold=Le;Le=Lb; Lb=Bhold;      
    Bhold=De; De=Db; Db=Bhold;
    Bhold=Se; Se=Sb; Sb=Bhold;
    Bhold=d; d=D; D=Bhold;
    %nib : need equation for band gap narrowing of p-type GaAs nib is Eg
    %narrowing effects for n-type GaAs
elseif pol=='pin'
    ZZZ=1; %no changes made to variables in MaterialParams.m
else
    disp(sprintf('Error: pin or nip polarity'))
end


%Absorption Coefficients
An=Aryanadj2; %Where do these come from - MaterialParams?????????????
Ap=Aryanadj2;
Ai=Aryanadj2;
%subscript w=window (a in Hovel) e=emitter (g in Hovel), b=base (p Hovel)
%absorption coefficients taken from Aguinaldo data p=pGaAs, n=nGaAs, and
%measured data (Slocum)
Q=B*Lw; %B is absorption in window layer
Y=-D/Lw;
Z=-Sw*tw/Lw;

S=Ap.*Le;
T=-d/Le;
U=-Se*te/Le;

E=Ap.*Lb;
Fa=Tb/Lb;
G=Sb*tb/Lb;

QEwAlgora=-(1-ref).*exp(-B.*D).*(Q./(Q.^2-1).*(Q-1+( (Z-Q).*exp(-Q.*Y)+(1-Z).*exp(-Y))./(cosh(Y)+Z*sinh(Y))));
QEeAlgora=-(1-ref).*exp(-B.*D-Ap.*d).*(S./(S.^2-1).*(S-1+( (U-S).*exp(-S.*T)+(1-U).*exp(-T))./(cosh(T)+U*sinh(T))));
QEbAlgora=(1-ref).*exp(-B.*D-Ap.*d-Ap.*Wd-Ap*wi).*(E./(E.^2-1).*(E-1+( (G-E).*exp(-E.*Fa)+(1-G).*exp(-Fa))./(cosh(Fa)+G*sinh(Fa))));

QEw=(1-ref).*B.*Lw./(B.^2.*Lw.^2-1).*( (B.*Lw+Sw.*tw./Lw.*(1-exp(-B.*D).*cosh(D./Lw))-exp(-B.*D).*sinh(D./Lw))./(Sw.*tw./Lw.*sinh(D./Lw)+cosh(D./Lw))-B.*Lw.*exp(-B.*D) );
QEe=(1-ref).*exp(-B.*D).*Ap.*Le./(Ap.^2.*Le.^2-1).*((Ap.*Le+Se.*te./Le.*(1-exp(-Ap.*d).*cosh(d./Le))-exp(-Ap.*d).*sinh(d./Le))./(Se.*te./Le.*sinh(d./Le)+cosh(d./Le))-Ap.*Le.*exp(-Ap.*d))+QEw./(Se.*te./Le.*sinh(d./Le)+cosh(d./Le)) ;
QEb=(1-ref).*exp(-B.*D).*An.*Lb./(An.*Lb+1).*exp(-An.*d).*exp(-An.*Wd); %from Hovel
QEscr=(1-ref).*exp(-B.*D).*(1-exp(-Ai.*(Wd+wi))).*exp(-An.*(d-wp)); %from Hovel (approximates d>>wp) nearly the same as Nelson

%from Nelson:
NQEe=exp(-B*D).*(1-ref).*Ap.*Le./(Ap.^2*Le^2-1).*( ( Se*Le/De+Ap*Le-exp(-Ap*(d-wp)).*(Se*Le./De*cosh((d-wp)/Le)+sinh((d-wp)/Le)))./(Se*Le/De*sinh((d-wp)/Le)+cosh((d-wp)/Le))-Ap.*Le.*exp(-Ap*(d-wp)))+QEw./(Se.*te./Le.*sinh(d./Le)+cosh(d./Le)) ;
NQEscr=(1-ref).*exp(-B.*D).*(1-exp(-Ai.*(Wd+wi))).*exp(-An.*(d-wp)); %contains attenuation from emitter layer
NQEb=(1-ref).*exp(-B.*D).*exp(-Ai.*wi).*An*Lb./(An.^2*Lb^2-1).*exp(-Ap*d-An*wn).*(An*Lb-(Sb*Lb/Db*(cosh((Tb-wn)/Lb)-exp(-An*(Tb-wn)))+sinh((Tb-wn)/Lb)+An.*Lb.*exp(-An*(Tb-wn)))./(Sb*Lb/Db*sinh((Tb-wn)/Lb)+cosh((Tb-wn)/Lb)));

%totals
QEtAlgora=QEeAlgora+QEwAlgora+QEbAlgora+NQEscr;
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


Lambda2=circshift(Lambda,1);
Lambda2(1)=[];
dLambda=Lambda(2:end)-Lambda2;
%jphi=X*sum(F.*dLambda.*(NSR(2:end)));% Integrated Current Density
%jphi09A161=X*sum(F.*dLambda.*(SR09A161(2:end)));% Integrated Current Density
%jphi10A018C3=X*sum(F.*dLambda.*(SR10A018C3(2:end)));% Integrated Current Density
%jphiam0=X*sum(F.*dLambda.*(NSR(2:end)));% Integrated Current Density
% Integrated Current Density of individual current components
%jphe=X*sum(F.*dLambda.*(NSRe(2:end)));
%jphb=X*sum(F.*dLambda.*(NSRb(2:end)));
%jphscr=X*sum(F.*dLambda.*(NSRscr(2:end)));
%Li=.3;
%widthi=.3;
%deltan=te*jphi/(q*d);
%desired output parameters
%FF=[0 0;0 0]
%Eff=[0 0;0 0]

%loop changing size
% for(j=1:5)
%     L=Li/j;
%     width=widthi/j;
%     area=length*width
% for(i=1:5)
    ps;
    jphi;
    jph=jphi*(1-ps);
  

    
%integration of Baseline SR (input variable:SR09A161)
% baseLambda2=circshift(baseLambda,1);
% baseLambda2(1)=[];
% basedLambda=baseLambda(2:end)-baseLambda2;
% baseJph=sum(F(2:end,spect).*basedLambda.*QEt(2:end));

%dark currents      same as Nelson, Algora. Hovel uses long base approx

J01e=q.*ni.^2./Na.*De./Le.*(sinh((d-wp)./Le)+Se.*Le./De.*cosh((d-wp)./Le))./(Se.*Le./De.*sinh((d-wp)./Le)+cosh((d-wp)./Le));
J01b=q.*nib.^2./Nd.*Db./Lb.*(sinh((Tb-wn)./Lb)+Sb.*Lb./Db.*cosh((Tb-wn)./Lb))./(Sb.*Lb./Db.*sinh((Tb-wn)./Lb)+cosh((Tb-wn)./Lb));
J1e=J01e*(exp(V/Vt)-1);
J1b=J01b*(exp(V/Vt)-1);
%J01h=q*ni^2/Nd*Db/Lb*coth( (Tb-d+Wd+wi)/Lb) Hovel base contribution
%assuming long base
J01=J01e+J01b;

J02dep=pi/2*q*nib*(Wd+wi)*C.*sinh(V./(2*Vt))./((Vbi-V)/(2*Vt));
J02perim=q*nib*sL*Lper/area;
J2dep=J02dep;
J2perim=J02perim.*(exp(V/(2*Vt))-1);
J1=J01.*(exp(V./Vt)-1);
J2=J02dep+J2perim;
Jdark=J1+J2;
Jt=jph-Jdark;  %total current 
J2al=pi/2*q*ni*(Wd+wi)*C/((Vbi-V)/(2*Vt))*(exp(V/(2*Vt))-1)+J02perim.*(exp(V/(2*Vt))-1);
Jtal=jph-J1-J2al;
Jtpos=(Jt+abs(Jt))/2;
%plot(V(1:112),J02(1:112))
%plot(V,Jtpos);
%semilogy(V,[J1 J2])
%max power, jsc, voc
power=Jt.*V; %[watts/cm2]      

Pmax=max(power); %[Watts/cm2]

n=1;
vmax=1;
jmax=1;
while( power(n)<Pmax )
    vmax=V(n);
    jmax=Jt(n);
    n=n+1;
end

k=1;
voc=1;

while(power(k)>=0)
   voc=V(k);
   k=k+1;
end

jsc=jph;

FFi=Pmax/(jsc*voc); %FF/Eff taking into account only grid shadowing
Effi=Pmax/(X*Psun); 

%F has 3 columns for AM0, 1.5 d and g
%make sure SR has correct units./formula
%put F in units of W./cm.^2./nm

%Grid stuff%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%power loss iterations

%intitial finger spacing guess
S1=(6*wf*vmax/(Rshe*jmax))^(1/3)*5;
numf=fix(L/S1);
S1=L/numf;
%intitial power losses
pe=Rshe/12*jmax/vmax*S1^2; %loss in emitter
prf=1/3*L^2*Rshm*jmax/vmax*S1/wf; %loss in fingers
pc=rhoc*jmax/vmax*S1/wf;
psf=wf/S1;
%iterate to get final S
tolerance=S1/100;
Si=S1;
sp=Si+2*tolerance;

c=0;

while(abs(sp-Si)>tolerance)  %loop iterates to minimize power loss, matches sheiman program
    Si=sp;
    sp=Si*2*(pe+psf)/(prf+4*pe+pc+psf);
    numf=fix(L/sp); %forces finger spacing to produce integer number of fingers
    sp=L/numf;
    %recalculate power loss..
    pe=Rshe/12*jmax/vmax*sp^2; 
    prf=1/3*width^2*Rshm*jmax/vmax*sp/wf; 
    pc=rhoc*jmax/vmax*sp/wf;    
    psf=wf/sp;
    c=c+1;
    numf;
    Matrix99(c,:)=[c numf sp pe prf pc psf];
    pt=pe+prf+pc+psf;
    abs(sp-Si);
%     if c>1000
%         break
%     end
end

%busbar
wb=L*width*sqrt(1/3*Rshm*jmax/vmax);
psb=wb/width;
prb=1/3*L^2*width/wb*Rshm*jmax/vmax;

%base
prbase=rhob*Tb/(L*width)*jmax/vmax;
%pri=rhoi*wi/(length*width)*jmax/vmax; %intrinsic region power loss..
prsub=rhosub*Tsub/(areaw)*jmax/vmax;
pbc=rhoc*jmax/vmax;  %back contact 
%total shadowing
ps=psf+psb;
%total resistive

pr=prf+prb+pe+pc+prbase+prsub+pbc;
rtotal=pr*vmax/jmax; %[ohm*cm^2]
Rtotal=rtotal/(area^2); %[ohm]
%num fingers

%end
% end
%add contact pad:
%resulting params

%area=length*width;
FF=FFi*(1-pr);
Pmax=Pmax*(1-pr);
Eff=Pmax/(X*Psun);
X;
isc=jsc*(1-ps); %is isc as modeled including fraction of reflected light from grid
%efficiency(j,1)=X
%efficiency(j,2)=Eff
Matrix=[sp,prb,prf,pe,prbase+prsub,psf,psb,L, width, area,jmax,vmax,jsc,voc,FF,Eff];
%pdark=1/jmax*(jmax-J01.*exp(vmax./Vt)-pi/2*q*ni*Wd*C.*2.*sinh(vmax./(2*Vt))./(q.*(Vbi-vmax)/(kb.*T))+q*ni*sL*Lper/area*exp(vmax/(2*Vt)))
%end
% re=Rshe/12*sp^2; 
% prf=1/3*width^2*Rshm*sp/wf; 
% pc=rhoc*sp/wf;    
% prb=1/3*L^2*width/wb*Rshm;
% prbase=rhob*Tb/(L*width);
% prsub=rhosub*Tsub/(areaw);
%  prbase+prsub

%series resistance effect on dark current
J1etest=J01e*(exp((V+rtotal*jmax)/Vt)-1);
J1btest=J01b*(exp((V+rtotal*jmax)/Vt)-1);
J02deptest=pi/2*q*nib*(Wd+wi)*C.*sinh((V+rtotal*jmax)./(2*Vt))./((Vbi-(V+rtotal*jmax))/(2*Vt));
J2perimtest=J02perim.*(exp((V+rtotal*jmax)/(2*Vt))-1);
Jdarktest=J1etest+J1btest+J02deptest+J2perimtest;
 Matrix2=[d,Na,Eff, FF, jsc, voc, pe,prf,prb,prbase+prsub,psf,psb];
 Matrix3=[Na,Nd,d, Tb, tf, wf, Eff, FF, jsc, voc, pe,prf,prb,prbase+prsub,psf,psb];
 NSRall=[NSR NSRb NSRe NSRscr];
 Matrix4=[Na tsrhe tsrhb tb te d jsc voc Pmax FF prf prb pe psf psb Eff  mue mub Le Lb];
 Matrix5=[Rshe pe Eff sp prf ps pc ];
 %Ploted output: SRw NSRe NSRb NSRscr NSR(From Nelson)
 output1=[Lambda QEw NQEe NQEscr NQEb SRw NSRe NSRb NSRscr NSR]; %Nelson SR outputs
%  output1=[Lambda QEwAlgora QEeAlgora QEscr QEbAlgora SRwAlgora SReAlgora SRbAlgora SRscr SRAlgora]; %Algora SR outputs
%  output1=[Lambda QEw QEe QEscr QEb SRw SRe SRb SRscr SR]; %Hovel SR outputs
 output2=[Effi FFi Eff FF Pmax jmax vmax voc jsc area J01e J01b J02dep(1) J02perim sp wb ps psb pc pbc pe prf prb prbase prsub numf]; %IV performance params and grid
 output3=[V Jt Jtpos J1e J1b J2dep J2perim power ]; %IV curve, dark current
 output=output1;
 switch(uu)
     case 'SR' 
         output=output1;
     case 'IV'
         output=output2;
     case 'dark' 
         output=output3;
 end
end
     