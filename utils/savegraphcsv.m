function savegraphcsv(G, filename, filedir)

if ~isdir(filedir)
    mkdir(filedir);
end
[a,b] = find(tril(G));
% Edges csv file
fid = fopen([filedir, filename],'w');
fprintf(fid, 'Source,Target,Type\n');
fprintf(fid, '%d,%d,Undirected\n', [a,b]');
fclose(fid);