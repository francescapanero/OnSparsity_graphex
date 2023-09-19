function out = omegafunc(x, y, pi, B)


% stochastic block-model graphon
% K = size(B, 1);
pi_sum = cumsum(pi);


[~,~,indx] = histcounts(x,[0,pi_sum]);
[~,~,indy] = histcounts(y,[0,pi_sum]);


% indx = find(x<pi_sum, 1);
% indy = find(y<pi_sum, 1);
ind = sub2ind(size(B), indx, indy);
out = B(ind);