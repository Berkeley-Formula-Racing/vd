function files = saveFigures(outDir,figs,opts)
% Writes open figures to disk, named from their figure Name.
%
%   saveFigures(dir)          every open figure, PNG at 200 dpi
%   saveFigures(dir,figs)     just these handles, in the order given
%   saveFigures(dir,[],opts)  every open figure, with options
%
% Pass an empty dir to do nothing and return {}. That is deliberate: the
% calling scripts expose the path as a setting, and '' is how you turn saving
% off without editing the call.
%
% Files are named <nn>_<figure Name>.<ext>. The number is the draw order, so
% the folder sorts the way the figures were created rather than
% alphabetically -- which for a sweep means the kernel and its overlay stay
% next to each other. It also keeps two figures that sanitise to the same
% name from overwriting one another.
%
% Existing files of the same name are overwritten. A sweep is expensive
% enough that losing its figures to a re-run is worth guarding against, so
% set opts.stamp = true to drop each run into its own timestamped subfolder.
%
% inputs  outDir: destination folder, created if missing. '' = do nothing.
%         figs:   figure handles, or [] for every open figure
%         opts:   .formats     cellstr, any of png jpg tiff pdf svg eps
%                              (default {'png'}). 'fig' also works and saves
%                              the reopenable MATLAB figure.
%                 .resolution  dpi for raster formats (default 200)
%                 .background  'white' (default) | 'none' | 'current' | RGB
%                 .prefix      prepended to every file name (default '')
%                 .number      draw-order prefix (default true)
%                 .stamp       timestamped subfolder (default false)
% output  files:  full paths written

if nargin < 1 || isempty(outDir) || strlength(string(outDir)) == 0
    files = {}; return
end
if nargin < 2, figs = []; end
if nargin < 3, opts = struct(); end

opts = withDefaults(opts,struct( ...
    'formats',{{'png'}}, 'resolution',200, 'background','white', ...
    'prefix','', 'number',true, 'stamp',false));
if ischar(opts.formats) || isstring(opts.formats)
    opts.formats = cellstr(opts.formats);
end

if isempty(figs)
    % findall rather than findobj: some figures are created with
    % HandleVisibility off and findobj would silently skip exactly those
    figs = findall(groot,'Type','figure');
    figs = sortByNumber(figs);
end
figs = figs(isgraphics(figs));
if isempty(figs)
    warning('saveFigures:noFigures','no open figures to save');
    files = {}; return
end

outDir = char(outDir);
if opts.stamp
    outDir = fullfile(outDir,char(datetime('now','Format','yyyy-MM-dd_HHmm')));
end
if ~isfolder(outDir)
    [ok,msg] = mkdir(outDir);
    if ~ok, error('saveFigures:mkdir','cannot create %s: %s',outDir,msg); end
end

files = {};
for i = 1:numel(figs)
    f = figs(i);
    stem = figName(f,i,opts);
    for k = 1:numel(opts.formats)
        ext  = lower(strtrim(opts.formats{k}));
        file = fullfile(outDir,[stem '.' ext]);
        if writeOne(f,file,ext,opts)
            files{end+1} = file; %#ok<AGROW>
        end
    end
end

fprintf('  saved %d file(s) to %s\n',numel(files),fullPath(outDir));
end


function ok = writeOne(f,file,ext,opts)
ok = true;
try
    switch ext
        case 'fig'
            savefig(f,file);
        case {'pdf','svg','eps'}
            exportgraphics(f,file,'ContentType','vector', ...
                'BackgroundColor',opts.background);
        otherwise
            exportgraphics(f,file,'Resolution',opts.resolution, ...
                'BackgroundColor',opts.background);
    end
catch err
    % Mostly for figures holding uicontrols -- the sim viewer and the track
    % builder are both like that. R2026a writes those and warns that the UI
    % components are dropped (use exportapp if you want them in the image);
    % older releases error outright. print goes through the ordinary figure
    % copy path either way, at the cost of honouring fewer options, so it is
    % the fallback rather than the default.
    try
        print(f,file,['-d' printDriver(ext)],sprintf('-r%d',opts.resolution));
    catch
        warning('saveFigures:failed','could not write %s: %s', ...
            file,err.message);
        ok = false;
    end
end
end


function d = printDriver(ext)
switch ext
    case 'jpg',  d = 'jpeg';
    case 'tiff', d = 'tiff';
    case 'eps',  d = 'epsc';
    otherwise,   d = ext;
end
end


function stem = figName(f,i,opts)
% Figure Name if it has one, otherwise the number MATLAB shows in the title
% bar. Everything that is not a letter, digit or underscore collapses to a
% single underscore, so 'ClA lap time attribution (autocross)' lands as
% ClA_lap_time_attribution_autocross and stays a legal file name on Windows.
nm = char(f.Name);
if isempty(strtrim(nm))
    if isempty(f.Number), nm = sprintf('figure_%d',i);
    else,                 nm = sprintf('figure_%d',f.Number);
    end
end
nm = regexprep(nm,'[^A-Za-z0-9]+','_');
nm = regexprep(nm,'^_+|_+$','');
if isempty(nm), nm = sprintf('figure_%d',i); end
if numel(nm) > 80, nm = nm(1:80); end

stem = [char(opts.prefix) nm];
if opts.number, stem = sprintf('%02d_%s',i,stem); end
end


function figs = sortByNumber(figs)
% findall returns newest first. Sort back into creation order so the file
% numbering matches the order the script drew them. Figures made with
% IntegerHandle off have no Number; they go last, in the order found.
n = arrayfun(@(h) numOrInf(h),figs);
[~,ord] = sort(n);
figs = figs(ord);
end


function n = numOrInf(h)
n = h.Number;
if isempty(n), n = inf; end
end


function s = withDefaults(s,d)
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(s,f{i}) || isempty(s.(f{i})), s.(f{i}) = d.(f{i}); end
end
end


function p = fullPath(p)
if ~contains(p,':') && ~startsWith(p,'\\'), p = fullfile(pwd,p); end
end
