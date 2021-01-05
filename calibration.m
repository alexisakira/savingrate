clear
clc;

set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultAxesTickLabelInterpreter','latex');
set(0,'DefaultLegendInterpreter', 'latex')
   
set(0,'DefaultTextFontSize', 14)
set(0,'DefaultAxesFontSize', 14)
set(0,'DefaultLineLineWidth',2)

colors = get(gca,'ColorOrder');

c1 = colors(1,:);
c2 = colors(2,:);

close all

%% use Saez & Zucman (2016) data to construct portfolio

year = xlsread('AppendixTables(Distributions).xlsx','TableB5b','A9:A108');
data001 = xlsread('AppendixTables(Distributions).xlsx','TableB5b','U9:Y108'); % top 0.01%
theta = mean((data001(:,1) + data001(:,4) + data001(:,5))./sum(data001,2));
% define risky portfolio by equity, business, and pension

%% estimate GARCH(1,1) model

GDP_Q = xlsread('RealPerCapitaGDP.xls','B12:B299');
g = mean(diff(log(GDP_Q)))/3; % monthly log GDP growth rate

% post 1947
Rfree = xlsread('PredictorData2018.xlsx','Monthly','K914:K1777');
infl = xlsread('PredictorData2018.xlsx','Monthly','L914:L1777');
CRSP_SPvw = xlsread('PredictorData2018.xlsx','Monthly','Q914:Q1777');
svar = xlsread('PredictorData2018.xlsx','Monthly','O914:O1777');

Rm = (1+CRSP_SPvw)./(1+infl); % real gross stock market return
Rex = (1+CRSP_SPvw)./(1+Rfree); % gross excess return
Rf = (1+Rfree)./(1+infl); % real gross risk-free rate

mu = mean(log(Rex)); % mean log excess return
sigma = std(log(Rex),1); % volatility

% discretize Gaussian distribution by Gauss-Hermite quadrature
J = 7; % number of states
[x, w] = GaussHermite(J);
PJ = w'/sqrt(pi); % probability
x = sqrt(2)*x';

tau_k = 0.25; % capital income tax rate
RX = 1 + (1-tau_k)*(exp(mu + sigma*x)-1)*theta; % realized asset return
Rf = exp(mean(log(Rf)) - g); % gross risk-free rate in detrended model

%% calibrate some parameters
gamma = 3; % relative risk aversion
eta = 1/25; % death rate (from length of one generation)
pd = 1- exp(-eta/12); % death probability (monthly)
tau_e = 0.4; % estate tax rate
zeta = 1.5; % US Pareto exponent

pUE = 1/3; % calibrate P(Unemployment -> Employment) by duration of unemployment
piU = 0.05; % unemployment rate
pEU = piU*pUE/(1-piU); % implied P(Employment -> Unemployment)
Pw = [1-pEU pEU; pUE 1-pUE]; % transition probability matrix for workers
yU = 0.2; % unemployment benefit (as fraction of wage)
ye = 2.5; % labor productivity of entrepreneurs (relative to workers)

phi = 0.115; % fraction of entrepreneurs ('active business owners' from Cagetti & De Nardi (2006))
pew = 1-exp(-0.025/12); % P(Entrepreneur -> Worker)
pwe = phi*pew/(1-phi); % P(Worker -> Entrepreneur)
PZ = [(1-pwe)*Pw pwe*ones(2,1); pew*[1-piU piU] 1-pew];
Z = size(PZ,1);

temp = [ones(2*Z,J); repmat(RX,Z,1)]; % gross excess returns conditional on survival
R = Rf*[temp (1-tau_e)*temp]; % matrix of returns
Y = repmat(repmat([1 yU ye]',1,2*J),Z,1); % matrix of incomes

PZt = @(t)([(1-phi*t/(1-phi))*Pw (phi*t/(1-phi))*ones(2,1); t*[1-piU piU] 1-t]);
func = @(t)((1-pd+pd*(1-tau_e)^zeta)*Rf^zeta*eigs(PZt(t).*repmat([1 1 dot(PJ,RX.^zeta)]',1,Z),1)-1);
t = fzero(func,[0,1]);
t = 1-exp(-0.02/12); % Gilchrist et.al report credit spread of 192 basis points
PZ = PZt(t);
