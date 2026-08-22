function name = setPlotFont(fig,requested)
% Applies the project plot font.
%
%   setPlotFont()      sets it session-wide on the graphics root, so EVERY
%                      figure created afterwards inherits it. setup_paths
%                      calls this, which is why plots come out in the project
%                      font without each plotting function asking.
%   setPlotFont(fig)   sets it on one figure only.
%
% Either way these are DEFAULT properties: they reach children created AFTER
% they are set, never existing ones. Per-figure callers must therefore call
% immediately after figure() and before any axes.
%
% MATLAB substitutes its own default when a FontName is not installed, and it
% does so silently: the figure renders, just in the wrong typeface. This
% resolves the name against listfonts first and warns once per session, so a
% missing font is visible instead of looking like a working one.
%
% 'Book' is a WEIGHT, not a family, and MATLAB picks weights through
% FontWeight rather than through the family name. Titles default to bold, so
% they are forced back to normal here -- otherwise a title asking for Futura
% Std Book renders in the bold face instead.
%
% inputs  fig:       figure handle, or omitted/[] for session-wide
%         requested: font name, or a cell array tried in order (optional).
%                    Defaults to the Futura Std Book family and its common
%                    alternative spellings.
% output  name:      the font actually applied

persistent warned

if nargin < 1 || isempty(fig), fig = groot; end
if nargin < 2 || isempty(requested)
    % Tried in order, first one that ACTUALLY RENDERS wins.
    %
    % "Futura Std Book" is listed by listfonts and by GDI+, is accepted by
    % FontName, reads straight back from the property -- and does not render.
    % Measured against a deliberately non-existent font name, it produces
    % identical metrics, i.e. MATLAB is substituting. Of the seven Futura Std
    % families Windows exposes, MATLAB can only render four: Condensed,
    % ExtraBold, Light and Medium. Book, Condensed ExtBd and Condensed Light
    % all fail. No amount of installing fixes it -- the family is installed
    % system-wide in C:\Windows\Fonts with an HKLM entry and still fails.
    %
    % "Futura Std" is not in listfonts at all, yet it renders -- as Medium,
    % confirmed pixel-identical. It is the nearest weight to Book that works
    % (Adobe's order is Light < Book < Medium), so it is the fallback.
    %
    % To get Book specifically, convert FuturaStdBook.otf to TrueType with
    % its own family name and install that; MATLAB renders it once it is a
    % family in its own right.
    requested = {'Futura Std Book','Futura Std','Futura Std Medium','Futura Std Light'};
end
if ischar(requested) || isstring(requested)
    requested = cellstr(requested);
end

% LISTFONTS IS NOT AUTHORITATIVE. It reports what Java AWT can enumerate,
% which is not the same list the figure renderer draws with. A font can be
% listed, be accepted by FontName, be reported straight back by the property
% -- and still be silently substituted at draw time. Measured on this
% machine: text set in "Futura Std Book" rendered pixel-identical to Times
% New Roman, while the property cheerfully read back "Futura Std Book".
%
% So the name is only a candidate. fontRenders() decides, by measuring text
% metrics against a font name that certainly does not exist: if the two match
% the renderer is substituting, whatever listfonts and the property claim.
% Walk the candidates and take the first that renders. resolve() alone is not
% enough -- it only proves a name is LISTED.
name = '';
for i = 1:numel(requested)
    if fontRenders(requested{i}), name = requested{i}; break, end
end
if ~isempty(name)
    if ~strcmpi(name,requested{1}) && isempty(warned)
        warning('setPlotFont:fellBack', ...
            ['"%s" does not render in MATLAB (it is listed, and silently ' ...
             'substituted). Using "%s" instead, the nearest candidate that ' ...
             'actually draws.\nTo get "%s" itself, convert its .otf to ' ...
             'TrueType under its own family name and install that.'], ...
             requested{1},name,requested{1});
        warned = true;
    end
    set(fig, 'DefaultAxesFontName',name, ...
             'DefaultTextFontName',name, ...
             'DefaultAxesTitleFontWeight','normal');
    return
end

name = resolve(listfonts,requested);

if isempty(name)
    % A font installed the usual way on Windows -- right click, Install --
    % goes to the PER-USER store under %LOCALAPPDATA%\Microsoft\Windows\Fonts
    % and registers under HKCU, which Java does not read. Registering it into
    % this JVM makes listfonts see it. Note this does NOT make the renderer
    % see it, which is exactly why the check below exists.
    registerUserFonts();
    name = resolve(listfonts,requested);
end

if ~isempty(name) && ~fontRenders(name)
    if isempty(warned)
        warning('setPlotFont:substituted', ...
            ['"%s" is listed but does NOT render -- MATLAB is substituting a ' ...
             'default and the figures will not be in it.\n' ...
             'This is the per-user font install: the figure renderer reads ' ...
             'C:\\Windows\\Fonts and HKLM, and a font installed just for you ' ...
             'lives in %%LOCALAPPDATA%%\\Microsoft\\Windows\\Fonts under HKCU. ' ...
             'No amount of code fixes it.\n' ...
             'FIX: right-click the .otf and choose "Install for all users" ' ...
             '(needs admin), then restart MATLAB.\n' ...
             'Falling back to %s so the figures are at least consistent.'], ...
             name,get(groot,'FactoryAxesFontName'));
        warned = true;
    end
    name = '';
end

if isempty(name)
    if isempty(warned)
        warning('setPlotFont:notInstalled', ...
            ['none of these fonts render: %s\nMATLAB substitutes its default ' ...
             'without saying so. If the font IS installed, it is probably ' ...
             'installed for you only -- re-install with right click > ' ...
             '"Install for all users". Otherwise check the FAMILY name, ' ...
             'which is often not the filename.'], strjoin(requested,', '));
        warned = true;
    end
    name = get(groot,'FactoryAxesFontName');
    return
end

% set only defaults that exist on every supported release; legends and
% colorbars inherit the axes font
set(fig, 'DefaultAxesFontName',name, ...
         'DefaultTextFontName',name, ...
         'DefaultAxesTitleFontWeight','normal');
end

function ok = fontRenders(name)
% Does this font actually reach the renderer, or is MATLAB substituting?
%
% Compares text metrics against a font name that cannot exist. Every
% unavailable name collapses onto the same substitute, so if the requested
% font measures identically to a made-up one, it is being substituted --
% regardless of what listfonts says or what the FontName property reads back.
%
% Metrics rather than pixels: an Extent query is milliseconds and needs no
% file, where rendering to PNG and diffing would cost a hundred times more
% for the same answer. Cached, so this runs once per font per session.
persistent cache
if isempty(cache), cache = struct(); end
key = matlab.lang.makeValidName(name);
if isfield(cache,key), ok = cache.(key); return, end

prev = get(groot,'CurrentFigure');
f = figure('Visible','off','Position',[10 10 200 100],'IntegerHandle','off', ...
    'HandleVisibility','off');
try
    ax = axes('Parent',f,'Visible','off');
    s  = 'Hamburgefonstiv 0123';
    tA = text(ax,0,0,s,'FontSize',30,'FontName',name);
    tB = text(ax,0,0,s,'FontSize',30,'FontName','ZZQQNoSuchFamily9');
    drawnow;
    eA = get(tA,'Extent'); eB = get(tB,'Extent');
    ok = max(abs(eA(3:4)-eB(3:4))) > 1e-6;
catch
    ok = true;   % a failed probe must not block the caller
end
if isgraphics(f), close(f); end
if ~isempty(prev) && isgraphics(prev), set(groot,'CurrentFigure',prev); end

cache.(key) = ok;
end

function name = resolve(avail,requested)
name = '';
for i = 1:numel(requested)
    j = find(strcmpi(avail,requested{i}),1);
    if ~isempty(j), name = avail{j}; return, end
end
end

function registerUserFonts()
% Loads every font in the Windows per-user font store into this JVM, so
% listfonts can see it. Registration is JVM-global and survives for the
% session, so this runs at most once.
persistent done
if ~isempty(done), return, end
done = true;

if ~ispc || ~usejava('jvm'), return, end
dirs = {fullfile(getenv('LOCALAPPDATA'),'Microsoft','Windows','Fonts')};

ge = java.awt.GraphicsEnvironment.getLocalGraphicsEnvironment();
for d = 1:numel(dirs)
    if ~isfolder(dirs{d}), continue, end
    files = [dir(fullfile(dirs{d},'*.otf')); dir(fullfile(dirs{d},'*.ttf'))];
    for k = 1:numel(files)
        p = fullfile(files(k).folder,files(k).name);
        try
            % TRUETYPE_FONT is correct for .otf too -- Java uses it for the
            % whole OpenType family, CFF outlines included
            jf = java.awt.Font.createFont(java.awt.Font.TRUETYPE_FONT,java.io.File(p));
            ge.registerFont(jf);
        catch
            % a font this JVM cannot parse is not worth failing a plot over
        end
    end
end
end
