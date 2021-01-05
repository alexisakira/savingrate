calibration

PJ = kron([1-pd pd],PJ); % redefine distribution of transitory shocks
%% define grid

temp = repmat(PJ,Z^2,1);
r = @(theta)(eigs(PZ.*(reshape(sum(temp.*(R.^theta),2),Z,Z)'),1));
r(1)

N = 100; % number of grid points
aMin = 1e-6; % minimum wealth
aMed = 10; % median grid point
aMax = 1e4; % maximum grid point
a=aMin;
b=aMax;
c=aMed;
% construct exponential grid on [a,b] with median point c
s = (c^2-a*b)/(a+b-2*c); % shift parameter
aGrid = exp(linspace(log(a+s),log(b+s),N))-s;
Nmed = floor(N/2);
aGrid(1:Nmed) = linspace(aMin,aMed,Nmed); % replace bottom half points by evenly spaced grid

MaxIter = 1e4;
tol = 1e-5;
delta = 0.04; % discount rate

gamma1 = 2; % relative risk aversion
gamma2 = 4;
beta1 = exp(-delta/12+(1-gamma1)*g); % discount factor
beta2 = exp(-delta/12+(1-gamma2)*g);

[ca1,cbar1] = getC(gamma1,PZ,PJ,beta1,R,Y,aGrid,MaxIter,tol);
[ca2,cbar2] = getC(gamma2,PZ,PJ,beta2,R,Y,aGrid,MaxIter,tol);

% verify Pareto exponent is correct
J = length(PJ);
G1 = R.*(1-kron(cbar1,ones(Z,J)));
G2 = R.*(1-kron(cbar2,ones(Z,J)));
zetaHat1 = getZeta(PZ,PJ,0,G1);
zetaHat2 = getZeta(PZ,PJ,0,G2);

% compute wealth distribution
gstjn1 = zeros(Z^2,N*J);
gstjn2 = zeros(Z^2,N*J);
for z = 1:Z
    for zhat = 1:Z
        zz = Z*(z-1)+zhat;
        for j = 1:J
            gstjn1(zz,N*(j-1)+1:N*j) = R(zz,j)*(aGrid-ca1(z,:)) + Y(zz,j);
            gstjn2(zz,N*(j-1)+1:N*j) = R(zz,j)*(aGrid-ca2(z,:)) + Y(zz,j);
        end
    end
end

[Q1,pi1] = getQ(PZ,PJ,0,1,aGrid,gstjn1,G1,zetaHat1); % compute stationary distribution
[Q2,pi2] = getQ(PZ,PJ,0,1,aGrid,gstjn2,G2,zetaHat2);
aDist1 = sum(reshape(pi1,N,Z),2);
aDist2 = sum(reshape(pi2,N,Z),2);

tailProb1 = flipud(cumsum(flipud(aDist1)));
tailProb2 = flipud(cumsum(flipud(aDist2)));

figure
loglog(aGrid,tailProb1,'--','Color',c1); hold on
loglog(aGrid,tailProb2,'-','Color',c1);
xlim([1,aMax])
legend(['$\gamma=$ ' num2str(gamma1)],['$\gamma=$ ' num2str(gamma2)])
xlabel('Asset')
ylabel('Tail probability')

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'fig_wDist','-dpdf')

% plot saving rates against wealth percentile

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
plot(1-tailProb1,sr1(3,:),'--','Color',c1); hold on
plot(1-tailProb2,sr2(3,:),'-','Color',c1);