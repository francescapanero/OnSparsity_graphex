function samples = texprnd(lambda, a, n, m)

% Samples from a left truncated exponential distribution
samples = -log( - rand(n, m).*exp(-lambda*a) + exp(-lambda*a) )./lambda;

end