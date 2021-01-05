calibration

PJ = kron([1-pd pd],PJ); % redefine distribution of transitory shocks
%% compute r(K(theta)) for theta = 0, 1, 1-gamma

temp = repmat(PJ,Z^2,1);
r = @(theta)(eigs(PZ.*(reshape(sum(temp.*(R.^theta),2),Z,Z)'),1));

Ngam = 101;
gamGrid = linspace(0,5,Ngam);
deltaGrid = 0*gamGrid;

delta0Grid = 12*((1-gamGrid)*g + log(r(0)));
delta1Grid = 12*((1-gamGrid)*g + log(r(1)));
for n=1:Ngam
    deltaGrid(n) = 12*((1-gamGrid(n))*g + log(r(1-gamGrid(n))));
end

figure
plot(gamGrid,delta0Grid,'--','Color',c1); hold on
plot(gamGrid,delta1Grid,':','Color',c1);
plot(gamGrid,deltaGrid,'Color',c2);
text(1,0.05,'$\bar{c}(z)>0$','HorizontalAlignment','center')
text(1,0.065,'$r(K(1-\gamma))<1$','HorizontalAlignment','center')
text(3,-0.015,'$\bar{c}(z)=0$','HorizontalAlignment','center')
text(3,0,'$r(K(1-\gamma))\ge 1$','HorizontalAlignment','center')
xlabel('$\gamma$ (risk aversion)')
ylabel('$\delta$ (annual discount rate)')
ylim([-0.1,0.1])
legend('$r(K(0))=1$','$r(K(1))=1$','$r(K(1-\gamma))=1$','Location','SW')

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'fig_rK','-dpdf')

%% example of consumption functions

N = 1000; % number of grid points
aMin = 1e-6; % minimum wealth
aMed = 10; % median grid point
aMax = 1e10; % maximum grid point
a=aMin;
b=aMax;
c=aMed;
% construct exponential grid on [a,b] with median point c
s = (c^2-a*b)/(a+b-2*c); % shift parameter
aGrid = exp(linspace(log(a+s),log(b+s),N))-s;
%Nmed = floor(N/2);
%aGrid(1:Nmed) = linspace(aMin,aMed,Nmed); % replace bottom half points by evenly spaced grid

MaxIter = 1e4;
tol = 1e-5;
delta = 0.04; % discount rate

%% low gamma
gamma1 = 2; % relative risk aversion
beta1 = exp(-delta/12+(1-gamma1)*g); % discount factor

[ca1,cbar1] = getC(gamma1,PZ,PJ,beta1,R,Y,aGrid,MaxIter,tol);

figure
plot(aGrid,ca1(1,:),'--','Color',c1); hold on
plot(aGrid,ca1(2,:),':','Color',c1);
plot(aGrid,ca1(3,:),'Color',c2);
xlim([0 20])
title(['$\gamma=$ ' num2str(gamma1)])
legend('Employed','Unemployed','Entrepreneur','Location','SE')
xlabel('Asset')
ylabel('Consumption')

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'fig_cf1_small','-dpdf')

figure
plot(aGrid,ca1(1,:),'--','Color',c1); hold on
plot(aGrid,ca1(2,:),':','Color',c1);
plot(aGrid,ca1(3,:),'Color',c2);
xlim([0,aMax])
title(['$\gamma=$ ' num2str(gamma1)])
legend('Employed','Unemployed','Entrepreneur','Location','SE')
xlabel('Asset')
ylabel('Consumption')

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'fig_cf1_large','-dpdf')

figure
loglog(aGrid,ca1(1,:)./aGrid,'--','Color',c1); hold on
loglog(aGrid,ca1(2,:)./aGrid,':','Color',c1);
loglog(aGrid,ca1(3,:)./aGrid,'Color',c2);
loglog(aGrid,cbar1*ones(1,N),'k:','LineWidth',1);
xlim([1e-2,aMax])
ylim([1e-5,1])
title(['$\gamma=$ ' num2str(gamma1)])
legend('Employed','Unemployed','Entrepreneur','Asymptotic MPC','Location','NE')
xlabel('Asset')
ylabel('Consumption rate')

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'fig_cr1','-dpdf')

%% high gamma

gamma2 = 4; % relative risk aversion
beta2 = exp(-delta/12+(1-gamma2)*g); % discount factor

[ca2,cbar2] = getC(gamma2,PZ,PJ,beta2,R,Y,aGrid,MaxIter,tol);

figure
plot(aGrid,ca2(1,:),'--','Color',c1); hold on
plot(aGrid,ca2(2,:),':','Color',c1);
plot(aGrid,ca2(3,:),'Color',c2);
xlim([0 20])
title(['$\gamma=$ ' num2str(gamma2)])
legend('Employed','Unemployed','Entrepreneur','Location','SE')
xlabel('Asset')
ylabel('Consumption')

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'fig_cf2_small','-dpdf')

figure
plot(aGrid,ca2(1,:),'--','Color',c1); hold on
plot(aGrid,ca2(2,:),':','Color',c1);
plot(aGrid,ca2(3,:),'Color',c2);
xlim([0,aMax])
title(['$\gamma=$ ' num2str(gamma2)])
legend('Employed','Unemployed','Entrepreneur','Location','SE')
xlabel('Asset')
ylabel('Consumption')

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'fig_cf2_large','-dpdf')

figure
loglog(aGrid,ca2(1,:)./aGrid,'--','Color',c1); hold on
loglog(aGrid,ca2(2,:)./aGrid,':','Color',c1);
loglog(aGrid,ca2(3,:)./aGrid,'Color',c2);
xlim([1e-2,aMax])
ylim([1e-5,1])
title(['$\gamma=$ ' num2str(gamma2)])
legend('Employed','Unemployed','Entrepreneur','Location','NE')
xlabel('Asset')
ylabel('Consumption rate')

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'fig_cr2','-dpdf')

% compute median saving rates
J = length(PJ);
jmed = ceil(J/4); % median index of j conditional on survival
sr1 = zeros(Z,N);
sr2 = zeros(Z,N);
for z = 1:Z
    zz = Z*(z-1) + z;
    Rplus = max(exp(g)*R(zz,jmed)-1,0); % net return excluding capital loss
    ahat1 = exp(g)*(R(zz,jmed)*(aGrid-ca1(z,:)) + Y(zz,jmed));
    sr1(z,:) = (ahat1 - aGrid)./(Rplus*(aGrid-ca1(z,:)) + exp(g)*Y(zz,jmed));
    ahat2 = exp(g)*(R(zz,jmed)*(aGrid-ca2(z,:)) + Y(zz,jmed));
    sr2(z,:) = (ahat2 - aGrid)./(Rplus*(aGrid-ca2(z,:)) + exp(g)*Y(zz,jmed));
end

figure
semilogx(aGrid,sr1(1,:),'--','Color',c1); hold on
semilogx(aGrid,sr1(2,:),':','Color',c1);
semilogx(aGrid,sr1(3,:),'-','Color',c2);
yline(0,'k:','LineWidth',1);
xlim([aMin,aMax])
ylim([-0.5,1])
title(['$\gamma=$ ' num2str(gamma1)])
legend('Employed','Unemployed','Entrepreneur','Location','SW')
xlabel('Asset')
ylabel('Saving rate')

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'fig_sr1','-dpdf')

figure
semilogx(aGrid,sr2(1,:),'--','Color',c1); hold on
semilogx(aGrid,sr2(2,:),':','Color',c1);
semilogx(aGrid,sr2(3,:),'-','Color',c2);
yline(0,'k:','LineWidth',1);
xlim([aMin,aMax])
ylim([-0.5,1])
title(['$\gamma=$ ' num2str(gamma2)])
legend('Employed','Unemployed','Entrepreneur','Location','SW')
xlabel('Asset')
ylabel('Saving rate')

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'fig_sr2','-dpdf')