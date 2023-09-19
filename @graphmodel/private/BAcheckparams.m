function BAcheckparams( n )
%BACHECKPARAMS checks the parameters of a Barabasi-Albert (BA) model

% Copyright (C) Francois Caron, University of Oxford and Adrien Todeschini, Inria
% caron@stats.ox.ac.uk
% adrien.todeschini@inria.fr
% September 2015
%--------------------------------------------------------------------------

if ~isnumeric(n) || (floor(n)-n)~=0
    error('Parameter n must be an integer');
end

end

