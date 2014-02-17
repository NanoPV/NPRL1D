function [X, Jsc, Voc, FFs, Efficiency]=concFn(F, dLambda, QE_t, Lambda, J0)

Xrow=logspace(0,3,100);
X=Xrow';

for i=1:length(X)
Jsc(i)=X(i)*sum(F.*dLambda.*QE_t.*Lambda./1240)*1000;
V=linspace(0,2,201);
J=J0*(exp(V./.026)-1);
Jlight=Jsc(i)-J*1000;
Voc(i)=.026*log(Jsc(i)/((J0)*1000));
voc=Voc/.026;
Rch=Voc/(Jsc/1000);
Rs=0.2;
rs=Rs/Rch;
pmax0=max(Jlight.*V);
FFo=pmax0/(Jsc(i)*Voc(i));
FFs(i)=FFo*(1-1.1*rs)+rs^2/5.4;%(voc-log(voc+.720))/(voc+1)%
pmaxs=FFs(i)*(Jsc(i)*Voc(i));
Efficiency(i)=pmaxs/(90*X(i))*100;
end

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