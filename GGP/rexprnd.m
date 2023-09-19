function out = rexprnd(lambda, a)

% REXPRND Sample from a right truncated distribution of pdf
% lambda * exp(-lambda*x) / (1-exp(-lambda*a)) si x in [0,a], 0 sinon

u = rand(size(lambda));
out = -1./lambda .*log(1 - u.*(1-exp(-lambda*a)));
end
