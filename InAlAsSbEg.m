

% Direct bandgap values for room temperature. From Vurgaftman.
EgInAs=0.4;
EgAlAs=3.0;   % 0 EgG EgIndir E1 E2
EgAlSb=2.3;
EgInSb=0.174;

% Ternary bowing. From Lumb.
BowInAlAs=0.7;
BowAlAsSb=0.8;
BowInAsSb=0.67;
BowInAlSb=0.43;

BowInAlAsSb=6.256;  % Quaternary bowing, from Lumb

%x=0.509; y=0.048;   % x% In, y% Sb
quatfrac=linspace(0,1,11);
x=0.5633.*(1-quatfrac); y=0.472.*quatfrac;
%x=0.515; y=0.040;
% Do all alloys at the same time:
%[x,y]=meshgrid(linspace(0,1,1000),linspace(0,1,1000));
quatfrac=linspace(0,1,11)
EgInAlAsSb=x.*(1-y)*EgInAs+(1-x).*(1-y)*EgAlAs+x.*y*EgInSb+(1-x).*y*EgAlSb-x.*(1-x).*(1-y).*BowInAlAs-x.*(1-x).*y*BowInAlSb-(1-x).*y.*(1-y)*BowAlAsSb-x.*y.*(1-y)*BowInAsSb-x.*y.*(1-x).*(1-y)*BowInAlAsSb

%imagesc(x(1,:),y(:,1),EgInAlAsSb);
hold on;
%contourvalues=0.1:0.2:max(max(EgInAlAsSb));
%[c,ch]=contour(x,y,EgInAlAsSb,contourvalues,'k');
%clabel(c,ch);

% Calculate and draw lattice-matched alloys
aInAs=6.0583;
aAlAs=5.6611;
aAlSb=6.1355;
aInSb=6.4794;
aInP=5.8697;

xLM_InAlAs=(aInP-aAlAs)/(aInAs-aAlAs);
yLM_AlAsSb=(aInP-aAlAs)/(aAlSb-aAlAs);

%plot([0 xLM_InAlAs],[yLM_AlAsSb 0],'k--');


set(gca,'YDir','normal');  % fix imagesc's flipping of the y axis

xlabel('% In');
ylabel('% Sb');