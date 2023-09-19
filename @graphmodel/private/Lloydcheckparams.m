function Lloydcheckparams( n, sig, c, d )
%LLOYDCHECKPARAMS checks the parameters of a Erdos-Renyi graph

% Copyright (C) Francois Caron, University of Oxford and Adrien Todeschini, Inria
% caron@stats.ox.ac.uk
% adrien.todeschini@inria.fr
% September 2015
%--------------------------------------------------------------------------

if ~isnumeric(n) || ~isnumeric(sig) || ~isnumeric(c) || ~isnumeric(d)
    error('Parameters must be numeric')
end
if n<=0 || sig<=0 || c<=0 || d<=0
    error('Parameters must be strictly positive');
end

end

