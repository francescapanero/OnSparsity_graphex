function [N, T] = GGPrnd(alpha, sigma, tau, T, maxiter)

%GGPrnd samples points of a generalized gamma process.
% [N, T] = GGPrnd(alpha, sigma, tau, T)
%
% Samples the points of the GGP with L�vy measure
%   alpha/Gamma(1-sigma) * w^(-1-sigma) * exp(-tau*w)
%
% For sigma>=0, it samples points above the threshold T>0 using the adaptive
% thinning strategy described in Favaro and Teh (2013).
% -------------------------------------------------------------------------
% INPUTS
%   - alpha: positive scalar
%   - sigma: real in (-Inf, 1)
%   - tau: positive scalar
%
% Optional inputs
%   - T: truncation threshold; positive scalar
%   - maxiter: maximum number of iterations for the adaptive thinning
%     strategy (defalut=1e8)
%
% OUTPUTS
%   - N: points of the GGP
%   - T: threshold
% -------------------------------------------------------------------------
% EXAMPLE
% alpha = 100; sigma = 0.5; tau = 1e-4;
% N = GGPrnd(alpha, sigma, tau);

% -------------------------------------------------------------------------
% Reference:
% S. Favaro and Y.W. Teh. MCMC for normalized random measure mixture
% models. Statistical Science, vol.28(3), pp.335-359, 2013.

% Copyright (C) Francois Caron, University of Oxford and Adrien Todeschini, Inria
% caron@stats.ox.ac.uk
% adrien.todeschini@inria.fr
% September 2015
%--------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%
%%% TODO: return expected number of connections not generated
%%% TODO: C++ with openmp implementation
%%%%%%%%%%%%%%%%%


% Check the parameters of the GGP
GGPcheckparams(alpha, sigma, tau);

%% Finite activity GGP
if sigma<-1e-8
    rate = exp( log(alpha) - log(-sigma) + sigma*log(tau) );
    K = poissrnd(rate);
    N = gamrnd(-sigma, 1/tau, K, 1);
    N = N(N>0);
    T = 0;
    return;
end

%% Infinite activity GGP
sigma = max(sigma, 0); %%% set sigma=0 if in [-1e-8, 0]

Njumps = [];
if nargin<4
    % set the threshold automatically so that we sample of the order Njumps jumps
    % Number of jumps of order alpha/sigma/Gamma(1-sigma) * T^{-sigma} for sigma>0
    % and alpha*log(T) for sigma=0
    if sigma>.1
        Njumps = 20000; % Expected number of jumps
        T = exp(1/sigma*(log(alpha) - log(sigma) - gammaln(1-sigma) - log(Njumps)));
    else
        T = 1e-10;
    end
elseif T<=0
    error('Threshold T must be strictly positive');
end
    
if isempty(Njumps)
    if sigma>1e-3
        Njumps = floor(exp(log(alpha) - log(sigma) -gammaln(1-sigma) - sigma*log(T)));
    else
        Njumps = floor(-alpha*log(T));
    end
end

if Njumps > 1e7
    warning('Expected number of jumps = %d - press key if you wish to continue', Njumps);
    pause
end

if nargin < 5
    maxiter = 1e8;
elseif maxiter<=0
    error('maxiter must be a strictly positive integer');
end

N = zeros(ceil(Njumps+3*sqrt(Njumps)), 1);
t = T;
k = 0;

if tau<1e-8
    %% case tau==0
    log_cst = log(alpha) - gammaln(1-sigma) - log(sigma);
    msigma = -sigma;
    msigmainv = -1/sigma;
    
    for k=1:maxiter
        log_r = log(-log(rand)) - log_cst;
        if log_r > msigma * log(t)
            completed = true;
            break;
        end
        t = exp( msigmainv * log( t^msigma - exp(log_r) ) );
        N(k) = t;
    end
else
    %% case tau>0
    % Adaptive thinning strategy
    log_cst = log(alpha)-gammaln(1-sigma)-log(tau);
    sigmap1 = 1+sigma;
    tauinv = 1/tau;
    
    for i=1:maxiter
        log_r = log(-log(rand)); % Sample exponential random variable e of unit rate
        log_G = log_cst-sigmap1*log(t)-tau*t;
        if log_r > log_G
            completed = true;
            break;
        end
        t_new = t-tauinv*log(1-exp(log_r-log_G));
        if log(rand) < sigmap1*(log(t)-log(t_new))
            k = k+1;
            N(k) = t_new;
        end
        t = t_new;
    end
end

N = N(1:k);

% If too many computions, we increase the threshold T and rerun
if ~completed
    T = T*10;
    warning('T too small - Its was value increased at %f', T);
    N = GGPrnd(alpha, sigma, tau, T);
end

end
