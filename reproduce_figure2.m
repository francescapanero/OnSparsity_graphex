%% Reproducing results of <https://doi.org/10.1017/apr.2022.75 On sparsity, power-law, and clustering properties of graphex processes>
% 
% 
% This Matlab script shows how to produce figure 2 of the paper <https://doi.org/10.1017/apr.2022.75 
% On sparsity, power-law, and clustering properties of graphex processes>. Data 
% are sampled from the generalized gamma process graph model of <https://doi.org/10.1111/rssb.12233 
% Sparse graphs using exchangeable random measures.>
% 
% Authors: <http://www.stats.ox.ac.uk/~caron/ François Caron> (University of 
% Oxford) and <https://francescapanero.github.io Francesca Panero> (London School 
% of Economics and Political Science).
% 
% Citation: Caron, F., Panero, F., & Rousseau, J. (2022). On sparsity, power-law, 
% and clustering properties of graphex processes. _Advances in Applied Probability_, 
% 1-43.
% 
% Tested on Matlab R2023a.
%% General settings

close all
clear all

% Add paths
addpath('./GGP/', './utils/');

% Set the fontsize
set(0,'DefaultAxesFontSize',14)

% where to save plots
saveplots = true;
saveworkspace = false;
if saveplots
    rep = './results/';
    if ~isdir(rep)
        mkdir(rep);
    end
end

% Set the seed
% rng('default')
%% Figure 2a: empirical degree distribution against theoretical asymptotic degree distribution (power-law)

figure
tau = 2;
sigma_all = [.2];
alpha = 1000;
for i=1:length(sigma_all)
    sigma = sigma_all(i);

    % sample graph from GGP model by Caron, Fox 2017
    obj = graphmodel('GGP', alpha, sigma, tau); 
    G = graphrnd(obj,1e-6); 

    deg = sum(G, 2);
    Nj = hist(deg, 1:max(deg));
    h = loglog(1:max(deg), Nj/size(G, 1), 'o');
    set(h, 'markersize', 4, 'color',  [.8, .3, .3],  'markerfacecolor', [.8, .3, .3])
    out = exp(log(sigma)+gammaln(h.XData-sigma) - gammaln(h.XData+1) - gammaln(1-sigma));
    hold on
    plot(h.XData, out, '--g', 'linewidth', 3, 'color',  [.3, .3, .8],  'markerfacecolor', [.8, .3, .3])
    legend('Empirical', 'Asymptotic')
    xlim([0,1000])
    ylim([1e-5, 1])
    xlabel('Degree $j$', 'interpreter', 'latex', 'fontsize', 16)
    ylabel('Degree distribution')
    box off
    ccdf_deg = cumsum(Nj, 'reverse');
end
if saveplots
    savefigs(gcf, 'powerlaw', rep);
end
%% Figure 2b: empirical clustering coefficients against theoretical asymptotic limits 

tau = 2; sigma = .2;

% Computing the asymptotic limits according to Proposition 11
rho = @(x) x.^(-1-sigma)/gamma(1-sigma).*exp(-tau*x);
fun1 = @(t) (1/sigma*((2*t+tau).^sigma - tau^sigma)).^2 .* rho(t);
fun2 = @(x,y,z) (1-exp(-2*x.*y)).*(1-exp(-2*x.*z)).*(1-exp(-2*y.*z)).* rho(x).*rho(y).*rho(z);

denom = integral(fun1, 0, Inf)
numer = integral3(fun2, 0, Inf, 0, Inf, 0, Inf);
global_lim = numer/denom

fun3 = @(y,z)  y.*z.*(1-exp(-2.*y.*z)).*rho(y).*rho(z);
numer = integral2(fun3, 0, Inf, 0, Inf);
m = tau^(sigma-1);
local_lim = numer/m^2

alpha_all =[20:20:80, 100:50:500];
nsamples = 10; % Number of graphs to sample
for j=1:nsamples
    % Sample limit graph
    obj = graphmodel('GGP', max(alpha_all), sigma, tau);
    [Gall, wtrue, wtrue_rem] = graphrnd(obj,1e-6);
    theta = max(alpha_all)*rand(size(wtrue));
    
    for i=1:length(alpha_all)
        alpha =alpha_all(i)

        G = Gall(theta<alpha, theta<alpha); % subset of graph corresponding to current alpha
        G = G-diag(diag(G));% remove self-edges
        deg = sum(G, 2); %Determine node degrees
        G = G(deg>1, deg>1);
        G = double(G);
        deg = sum(G, 2); %Determine node degrees
        cn = diag(G*triu(G)*G); %Number of triangles for each node
        %The local clustering coefficient of each node
        c = zeros(size(deg)); c = 2 * cn ./ (deg.*(deg - 1));
        c2 = zeros(length(deg), 1);
        for k=2:max(deg)
            ind = find(deg==k);
            c2(k) = mean(c(ind));
        end
        
        local_clust(i, j) = mean(c(deg>1));
        global_clust(i, j) = full(sum(cn)/sum(deg(deg>1).*(deg(deg>1)-1)/2));
        
    end
end

figure
plot(alpha_all, global_clust, 'r')
hold on
plot(alpha_all,global_lim*ones(size(alpha_all)), '--r', 'linewidth', 3)
plot(alpha_all, local_clust, 'b')
plot(alpha_all, local_lim*ones(size(alpha_all)), '--b', 'linewidth', 3)
xlabel('Graph size \alpha')
ylabel('Clustering coefficients')
ylim([0,0.3])
if saveplots
    savefigs(gcf, 'clustering', rep);
end
%% Figure 2c: local clustering coefficient per degree for large alpha, empirical vs asymptotic

alpha = 2000;

rho = @(x) x.^(-1-sigma)/gamma(1-sigma).*exp(-tau*x);
fun = @(x,y) (1-exp(-2.*x.*y)).*rho(x).*rho(y);
b2 = integral2(fun, 0, Inf, 0, Inf);

obj = graphmodel('GGP', alpha, sigma, tau);
G = graphrnd(obj,1e-6);

G = G-diag(diag(G));% remove self-edges
deg = sum(G, 2); %Determine node degrees
G = G(deg>1, deg>1);
G = double(G);
deg = sum(G, 2); %Determine node degrees
cn = diag(G*triu(G)*G); %Number of triangles for each node
%The local clustering coefficient of each node
c = zeros(size(deg)); c = 2 * cn ./ (deg.*(deg - 1));
c2 =zeros(length(deg), 1);
for k=2:max(deg)
    ind = find(deg==k);
    c2(k) = mean(c(ind));
end

figure
semilogx(c2,'.', 'color', [.8,.3,.3], 'markersize', 10)
hold on
semilogx(local_lim*ones(size(c2)), '--b', 'linewidth', 3, 'color', [.3,.3,.8])
xlim([2, length(c2)])
xlabel('Degree $j$', 'interpreter', 'latex')
ylabel('Local clustering coefficients')
legend('Empirical', 'Asymptotic')
if saveplots
    savefigs(gcf, 'localclusteringj', rep);
end