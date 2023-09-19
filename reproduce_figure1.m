% Reproducing Figure 1 in On sparsity, power-law, and clustering properties of graphex
% processes (https://doi.org/10.1017/apr.2022.75)
% Illustration of the paper to show the construction of the model from Kallenberg
% (See Caron, Fox (2017) "Sparse graphs using exchangeable random measures" https://doi.org/10.1111/rssb.12233 
% for a thorough explanation)

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
addpath('./utils/');

% colour
basecol = [139,0,0]/255;
nmap = 50;
map = [linspace(1,basecol(1), 50)', linspace(1,basecol(2), 50)', linspace(1,basecol(3), 50)'];

rng(13) 
alpha_par = 1;
Tmax = 6;
N = poissrnd(Tmax*alpha_par);
theta = alpha_par*rand(N, 1);
Theta = Tmax*rand(N,1);

% Figure 1a: point process representation
h = figure
plot(theta, Theta, 'o', 'markersize', 10, 'markerfacecolor', basecol, 'MarkerEdgeColor', 'none')
box off
hold on
xlim([0,1.2])
ylim([0,Tmax])
plot([1,1],[0,Tmax], '--k', 'linewidth', 2)
plot([0, theta(1)], [Theta(1), Theta(1)], '--k')
plot([0, theta(2)], [Theta(2), Theta(2)], '--k')

plot([theta(1), theta(1)], [0, Theta(1)], '--k')
plot([theta(2), theta(2)], [0, Theta(2)], '--k')
set(gca, 'fontsize', 20)
if theta(1)<theta(2)
    set(gca, 'XTick', [theta(1), theta(2), 1], 'XTickLabel', {'\theta_i','\theta_j', '\alpha'})
else
    set(gca, 'XTick', [theta(2), theta(1), 1], 'XTickLabel', {'\theta_j','\theta_i', '\alpha'})
end
if Theta(1)<Theta(2)
    set(gca, 'YTick', [Theta(1),Theta(2)], 'YTicklabel', {'\vartheta_i', '\vartheta_j'})
else
    set(gca, 'YTick', [Theta(1),Theta(2)], 'YTicklabel', {'\vartheta_j', '\vartheta_i'})
end
saveas(gca, 'results/pointprocess', 'png')

inv_levy = @(y,sigma) (y*sigma*gamma(1-sigma)).^(-1/sigma);
M = @(x,y,sigma) 1-exp(-2*inv_levy(x,sigma).*inv_levy(y,sigma));

% Figure 1b: graphon W
h = figure

[X,Y] = meshgrid(0:.05:Tmax,0:.05:Tmax);
Z = M(X,Y,.5);
contourf(X,Y,Z, 20,'EdgeColor','none')
set(gca, 'fontsize', 20)
alpha .1
grid off
colormap(map)
if Theta(1)<Theta(2)
    set(gca, 'XTick', [Theta(1),Theta(2)], 'XTickLabel', {'\vartheta_i','\vartheta_j'})
    set(gca, 'YTick', [Theta(1),Theta(2)], 'YTickLabel', {'\vartheta_i','\vartheta_j'})
else
    set(gca, 'XTick', [Theta(2),Theta(1)], 'XTickLabel', {'\vartheta_j','\vartheta_i'})
    set(gca, 'YTick', [Theta(2),Theta(1)], 'YTickLabel', {'\vartheta_j','\vartheta_i'})
end
box off
hold on
plot([0, Theta(1)], [Theta(2), Theta(2)], '--k')

plot([Theta(1), Theta(1)], [0, Theta(2)], '--k')
plot(Theta(1), Theta(2), 'ro', 'markersize', 10, 'markerfacecolor', 'b', 'markeredgecolor', 'none')
saveas(gca, 'results/W', 'png')
