function savefigs(h, filename, dirpath, format, siz)
% SAVEFIGS save matlab figures
% 
% Optional Inputs:
% - h: figure handle vector
% - filename: output filename prefix (adds numbers suffix if sevral figures)
% - dirpath: output directory
% - format: string or cell of strings. file formats
% - siz: vector [width height] in pixels

if nargin<1
    h = gcf;
end
if nargin<2
    filename = 'fig';
end
if nargin<3
    dirpath = '.';
end
if nargin<4
    format = {'png', 'epsc2', 'fig'};
end

if ~iscell(format)
    format = {format};
end

if ~isdir(dirpath)
    mkdir(dirpath);
end

suffix = '';

for i=1:numel(h)
    if numel(h)>1
        suffix = sprintf('_%02d', i);
    end
    
    if nargin>=5
        pos = get(h(i), 'Position');
        pos(3:4) = siz;
        set(h(i), 'Position', pos)
        set(h(i), 'PaperUnit', 'points')
        set(h(i), 'PaperPosition', [1,1,siz])
        set(h(i), 'PaperSize', siz)
    end

    for j=1:numel(format)
        saveas(h(i), fullfile(dirpath, [filename, suffix]), format{j})
    end
end
