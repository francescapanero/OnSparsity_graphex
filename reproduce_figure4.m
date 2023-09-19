% Reproducing Figure 4 in On sparsity, power-law, and clustering properties of graphex
% processes (https://doi.org/10.1017/apr.2022.75)

% Authors: <http://www.stats.ox.ac.uk/~caron/ François Caron> (University of 
% Oxford) and <https://francescapanero.github.io Francesca Panero> (London School 
% of Economics and Political Science).
% 
% Citation: Caron, F., Panero, F., & Rousseau, J. (2022). On sparsity, power-law, 
% and clustering properties of graphex processes. _Advances in Applied Probability_, 
% 1-43.
% 
% Tested on Matlab R2023a.


close all
clear all
set(0,'DefaultAxesFontSize',14)
addpath('./utils/');

% set color
basecol = [139,0,0]/255;
nmap = 50;
map =[linspace(1,basecol(1), 50)', linspace(1,basecol(2), 50)', linspace(1,basecol(3), 50)'];

% set seed
rng(0);

% Graphon defining the sparse stochastic blochmodel (see Example 2 in the paper)
stepsize = .01;
[x, y] = meshgrid(0:stepsize:1, 0:stepsize:1);
pi = [.5, .3, .2];
B = [.7, .1, .1;
    .1, .5, .05;
    .1, .05, .9];
z = omegafunc(x, y, pi, B);

figure
imagesc(0:stepsize:1, 0:stepsize:1, z)
set(gca,'YDir','normal')
colormap(map)
xlabel('x', 'fontsize', 16)
ylabel('y', 'fontsize', 16)
caxis([0,1])
% Figure 4a
saveas(gca, 'results/graphon.png', 'png')
saveas(gca, 'results/graphon.eps', 'epsc2')

sigma = .8;

maxval = 6;
stepsize=.01;
[x, y] =meshgrid(0:stepsize:maxval, 0:stepsize:maxval);
z = eta(x, y, sigma);

figure
imagesc(0:stepsize:maxval, 0:stepsize:maxval, z)
set(gca,'YDir','normal')
colormap(map)
xlabel('x', 'fontsize', 16)
ylabel('y', 'fontsize', 16)
caxis([0,1])
% Figure 4b
saveas(gca, 'results/sparsegraphon.png', 'png')
saveas(gca, 'results/sparsegraphon.eps', 'epsc2')


W = @(u1, v1, u2, v2) omegafunc(v1, v2, pi, B).*eta(u1, u2, sigma);

alpha = 50;
trunc = 100; % truncation for vartheta
K = poissrnd(trunc*alpha);
vartheta = trunc*rand(K, 1);
v = rand(K, 1);
z = zeros(K);
for i=1:K-1
        z(i, i+1:K) = rand( K-i, 1)<W(vartheta(i)*ones(K-i,1), v(i)*ones(K-i,1), vartheta(i+1:K), v(i+1:K));
end
z = z + z' ;
ind = find(sum(z));
z = z(ind, ind);
vartheta = vartheta(ind);
size(z, 1)

% save graph as csv
% savegraphcsv(z, 'results/sparseSBM.csv', './');    

% Figure 4d: empirical degree distribution
figure
h = plot_degree(z, 'o')
set(h, 'markersize', 6, 'color',  basecol,  'markerfacecolor', basecol)
saveas(gca, 'results/degreedistribution.eps', 'epsc2' )
saveas(gca, 'results/degreedistribution.png', 'png' )

