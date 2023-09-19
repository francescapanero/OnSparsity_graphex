function print_features( filename, w, ind_features, featnames, meta, fnames, formats, n )

p = size(w,2);
if nargin<3
    ind_features = 1:p;
end
if nargin<4 || isempty(featnames)
    featnames = arrayfun(@(x) ['FEATURE ' num2str(x)], 1:p, 'uniformoutput', false);
end
if nargin<5
    meta = struct;
end
if nargin<6
    fnames = fieldnames(meta);
end
if nargin<7
    formats = repmat({'%s,'}, numel(fnames), 1);
end
if nargin<8
    n = 10;
end

fids = [1 fopen(filename, 'w')];
for k=1:p
    ind_k = ind_features(k);
    [~, ind] = maxn(w(:,ind_k), n);
    mfprintf(fids, '& %%--- %s ---\n', featnames{k})
    for j=1:n
%         mfprintf(fids, '%d. ', j)
        for i=1:numel(fnames)
            mfprintf(fids, [formats{i} ' '],meta.(fnames{i}){ind(j)})
        end
        mfprintf(fids, '\n')
    end
end
fclose(fids(2));
end


function [C, I] = maxn(A, n)
C = zeros(n,1);
I = zeros(n,1);
for i=1:n
    [C(i), I(i)] = max(A);
    A(I(i)) = -inf;
end
end
