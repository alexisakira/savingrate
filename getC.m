function [ca,cbar] = getC(gamma,PZ,PJ,beta,R,Y,aGrid,MaxIter,tol,ca0)
% compute consumption function by policy function iteration
% speed up using variant of endogenous grid point method

% gamma: relative risk aversion
% PZ: transition probability matrix of Markov states
% PJ: probability matrix of transitionary states
% beta: matrix of discount factors
% R: matrix of asset returns
% Y: matrix of income
% aGrid: asset grid
% MaxIter: maximum number of iterations
% tol: error tolerance
% ca0: initial guess

%% some error checking
if gamma <= 0
    error('gamma must be positive')
end
if any(PZ < 0)
    error('P must be nonnegative')
end
Z = size(PZ,1); % number of Markov states
J = size(PJ,2); % number of transitionary states
if size(PJ,1) == 1
    PJ = repmat(PJ,Z^2,1); % if PJ row vector, then assume iid
end

if any(beta <= 0)
    error('beta must be positive')
end
if numel(beta) == 1
    beta = beta*ones(Z^2,J); % if beta scalar, then assume constant discount factor
end

if any(R <= 0)
    error('R must be positive')
end
if size(R,1) ~= Z^2
    error('R has incorrect rows')
end
if size(R,2) ~= J
    error('R has incorrect columns')
end

if any(Y <= 0)
    error('Y must be positive')
end
if size(Y,1) ~= Z^2
    error('Y has incorrect rows')
end
if size(Y,2) ~= J
    error('Y has incorrect columns')
end

if any(diff(aGrid) <= 0)
    error('aGrid must be strictly increasing')
end
if aGrid(1) <= 0
    error('elements of aGrid must be positive')
end
if size(aGrid,1) > size(aGrid,2)
    aGrid = aGrid'; % convert to row vector
end
N = length(aGrid); % number of grid points
% set default if not provided
if nargin < 8
    MaxIter = 1e4;
end
if nargin < 9
    tol = 1e-5;
end

%% check existence of solution and compute asymptotic MPC
% define matrix K(theta)
K = @(theta)(PZ.*(reshape(sum(PJ.*beta.*(R.^theta),2),Z,Z)'));
K0 = K(0);
K1 = K(1);
K = K(1-gamma);

% spectral radius
r0 = eigs(K0,1);
r1 = eigs(K1,1);
r = eigs(K,1);

fprintf('r0 = %0.10f\n',r0)
fprintf('r1 = %0.10f\n',r1)
fprintf('r = %0.10f\n',r)

% check existence condition
if (r0 >= 1)||(r1 >= 1)
    disp('Nonexistence')
    return
end

% compute asymptotic MPC
if r >= 1
    cbar = zeros(Z,1);
else
    F =@(x)(1 + (K*x).^(1/gamma)).^gamma; % fixed point mapping
    m0 = (1 - r^(1/gamma))*ones(Z,1); % initial guess
    g =@(m)((m - F(m.^(-gamma)).^(-1/gamma))/norm(m0)); % convert fixed point equation in space of MPC
    options = optimoptions('fsolve','FunctionTolerance',1e-10,'OptimalityTolerance',1e-10,'Display','off');
    cbar = fsolve(g,m0,options);
end

%% policy function iteration

EY = sum(PZ.*(reshape(sum(PJ.*Y,2),Z,Z)'),2); % expected income

ca_bind = repmat(aGrid,Z,1); % consumption function with binding borrowing constraint
if nargin == 10 % initial guess provided
    ca = ca0;
else
    ca = EY + max(cbar,1e-4)*aGrid; % initialize consumption function using asymptotic MPC
    ca = min(ca,ca_bind); % impose borrowing constraint
end

% policy function iteration
for i = 1:MaxIter
    ca_new = zeros(Z,N);
    for z = 1:Z % iterate over current state
        expect = 0*aGrid; % conditional expectation
        for zhat = 1:Z % iterate over next state
            for j = 1:J % iterate over transitory state
                phat = PZ(z,zhat)*PJ(Z*(z-1)+zhat,j); % conditional probability
                bhat = beta(Z*(z-1)+zhat,j); % discount factor
                Rhat = R(Z*(z-1)+zhat,j); % return
                Yhat = Y(Z*(z-1)+zhat,j); % income
                ahat = Rhat*(aGrid - ca(z,:)) + Yhat; % next period asset
                ahat = cummax(ahat); % force it to be an increasing function
                chat = interp1([0 aGrid],[0 ca(zhat,:)],ahat,'linear','extrap'); % linear interpolation
                chat = max(chat,1e-10); % force it to be nonnegative
                expect = expect + phat*bhat*Rhat*chat.^(-gamma); % update expectation
            end
        end
        ca_new(z,:) = expect.^(-1/gamma); % update consumption using Euler equation
    end
    ca_new = (ca + ca_new)/2; % don't update all the way (for numerical stability)
    ca_new = min(ca_new,ca_bind); % impose borrowing constraint
    %ca_new = cummax(ca_new); % force it to be an increasing function
    error_i = max(max(abs(ca_new./ca - 1))); % relative error in last step
    if error_i < tol
        disp('Converged!')
        break
    else
        ca = ca_new;
    end
end

fprintf('Number of iteration = %0.0f, error = %0.10f\n',i,error_i)

end

