function ERcheckparams( n, p )
%ERCHECKPARAMS checks the parameters of a Erdos-Renyi graph

% Copyright (C) Francois Caron, University of Oxford and Adrien Todeschini, Inria
% caron@stats.ox.ac.uk
% adrien.todeschini@inria.fr
% September 2015
%--------------------------------------------------------------------------

if ~isnumeric(n) || (floor(n)-n)~=0
    error('First parameter n must be an integer');
end
if ~isnumeric(p) || p>1 || p<0
    error('Second parameter p must be a real in (0,1)')
end

end

