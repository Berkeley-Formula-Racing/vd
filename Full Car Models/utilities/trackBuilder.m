function trackBuilder(imgFile)
% Interactive track builder: trace a course over a photo or map and write the
% .mat file the lap sim reads. Replaces drawing a spline in CAD and exporting
% by hand.
%
%   trackBuilder                % opens a file picker for the image
%   trackBuilder('michigan.png')
%
% WHAT THE LAP SIM ACTUALLY WANTS
% Events2 loads a .mat holding exactly two row vectors:
%   arclength  cumulative distance along the centreline, metres, increasing
%   curvature  signed 1/radius at each of those stations, 1/metres
% Nothing else. There are no x,y coordinates in the file -- the sim only ever
% needs how far you have gone and how tight it is where you are. Sign
% distinguishes left from right; the solver takes the magnitude when it turns
% curvature into a corner radius.
%
% WHY THIS IS NOT JUST "FIT A SPLINE"
% Curvature is the SECOND derivative of position, so it amplifies noise
% savagely. Clicking points a couple of pixels off a smooth line produces
% curvature that wobbles across zero every few centimetres, and curvature_apexes
% calls findpeaks with no prominence filter -- so every one of those wobbles
% becomes an "apex" the solver caps speed at. Measured on the existing files:
%
%   michigantrack2024   1615 apexes over 703 m  (one every 0.4 m)
%                       1547 of them at R > 200 m, i.e. on straights
%   2024endurancetrack   663 apexes over 2061 m (one every 3.1 m)
%
% Those are not corners. This is why the smoothing control exists, and why
% the readout shows the apex count the SOLVER will find rather than some
% internal number: tune until the apex count matches the number of real
% corners you can see on the image.
%
% HOW TO USE
%   1  Load image        the track photo, map screenshot, or course diagram
%   2  Set scale         click two points a known distance apart, type the
%                        distance in metres. Nothing works before this.
%   3  Add points        click along the CENTRELINE, Enter when done. You can
%                        add more later; they append to the end.
%   4  Closed loop       tick for endurance/autocross loops that rejoin. The
%                        spline is wrapped so the seam is smooth rather than
%                        showing up as a fake corner.
%   5  Smoothing         drag until the apex count settles near the number of
%                        real corners. Watch the curvature plot: real corners
%                        are broad humps, noise is grass.
%   6  Save track        writes arclength + curvature to a .mat.
%
% Load track re-opens an existing .mat to inspect it (curvature and apex
% count only -- the file has no geometry to draw).

if nargin < 1, imgFile = ''; end

S = struct();
S.img      = [];
S.imgFile  = '';
S.scale    = [];        % metres per pixel
S.pts      = zeros(0,2);% clicked points, pixels
S.closed   = false;
S.smooth   = 0.5;       % 0 = heavy smoothing, 1 = interpolate the clicks
S.nOut     = 100000;    % output stations, matching the existing files
S.fit      = [];

S.fig = figure('Name','trackBuilder','NumberTitle','off', ...
    'Position',[80 80 1280 720],'Color',[0.94 0.94 0.94]);
setPlotFont(S.fig);

S.axImg = axes('Parent',S.fig,'Units','normalized','Position',[0.04 0.30 0.52 0.64]);
title(S.axImg,'load an image to start'); axis(S.axImg,'off');
S.axCur = axes('Parent',S.fig,'Units','normalized','Position',[0.63 0.55 0.34 0.38]);
xlabel(S.axCur,'distance (m)'); ylabel(S.axCur,'curvature (1/m)'); grid(S.axCur,'on');
S.axRad = axes('Parent',S.fig,'Units','normalized','Position',[0.63 0.30 0.34 0.16]);
xlabel(S.axRad,'distance (m)'); ylabel(S.axRad,'radius (m)'); grid(S.axRad,'on');

b = @(x,str,cb) uicontrol(S.fig,'Style','pushbutton','String',str, ...
    'Units','normalized','Position',[x 0.16 0.115 0.055],'Callback',cb);
b(0.04,'Load image',   @(~,~)cbLoadImage(S.fig));
b(0.17,'Set scale',    @(~,~)cbScale(S.fig));
b(0.30,'Add points',   @(~,~)cbAdd(S.fig));
b(0.43,'Edit (drag)',  @(~,~)cbEdit(S.fig));
b(0.56,'Undo point',   @(~,~)cbUndo(S.fig));
b(0.69,'Save track',   @(~,~)cbSave(S.fig));
b(0.82,'Load track',   @(~,~)cbLoadTrack(S.fig));

S.hClosed = uicontrol(S.fig,'Style','checkbox','String','closed loop', ...
    'Units','normalized','Position',[0.04 0.10 0.12 0.04], ...
    'BackgroundColor',[0.94 0.94 0.94],'Callback',@(~,~)cbClosed(S.fig));

uicontrol(S.fig,'Style','text','String','smoothing', ...
    'Units','normalized','Position',[0.17 0.10 0.08 0.035], ...
    'BackgroundColor',[0.94 0.94 0.94],'HorizontalAlignment','left');
S.hSmooth = uicontrol(S.fig,'Style','slider','Min',0,'Max',1,'Value',S.smooth, ...
    'Units','normalized','Position',[0.25 0.105 0.30 0.03], ...
    'Callback',@(h,~)cbSmooth(S.fig,h));

S.hStat = uicontrol(S.fig,'Style','text','String','', ...
    'Units','normalized','Position',[0.04 0.015 0.93 0.075], ...
    'BackgroundColor',[0.94 0.94 0.94],'HorizontalAlignment','left','FontSize',9);

guidata(S.fig,S);
if ~isempty(imgFile), loadImage(S.fig,imgFile); end
status(S.fig,'Load an image, then Set scale before anything else.');
end

%% ------------------------------------------------------------------ actions
function cbLoadImage(fig)
[f,p] = uigetfile({'*.png;*.jpg;*.jpeg;*.tif;*.bmp','Images'},'Track image');
if isequal(f,0), return, end
loadImage(fig,fullfile(p,f));
end

function loadImage(fig,file)
S = guidata(fig);
try
    S.img = imread(file);
catch ME
    status(fig,sprintf('could not read image: %s',ME.message)); return
end
S.imgFile = file; S.pts = zeros(0,2); S.scale = []; S.fit = [];
guidata(fig,S); redraw(fig);
status(fig,sprintf('loaded %s (%d x %d px). Now Set scale.', ...
    file,size(S.img,2),size(S.img,1)));
end

function cbScale(fig)
S = guidata(fig);
if isempty(S.img), status(fig,'load an image first'); return, end
status(fig,'click TWO points a known distance apart ...');
axes(S.axImg);
[x,y] = ginput(2);
if numel(x) < 2, status(fig,'cancelled'); return, end
px = hypot(diff(x),diff(y));
a = inputdlg(sprintf('Distance between those points, in METRES\n(%.1f px apart):',px), ...
    'Set scale',1,{'10'});
if isempty(a), status(fig,'cancelled'); return, end
m = str2double(a{1});
if ~isfinite(m) || m <= 0, status(fig,'distance must be a positive number'); return, end
S.scale = m/px;
guidata(fig,S); refit(fig);
status(fig,sprintf('scale set: %.4f m/px (%.1f px = %.2f m). Now Add points.',S.scale,px,m));
end

function cbAdd(fig)
S = guidata(fig);
if isempty(S.img),   status(fig,'load an image first'); return, end
if isempty(S.scale), status(fig,'set the scale first'); return, end
status(fig,'click along the CENTRELINE; press Enter to finish');
axes(S.axImg);
[x,y] = ginput;
if isempty(x), status(fig,'no points added'); return, end
S.pts = [S.pts; x(:) y(:)];
guidata(fig,S); refit(fig);
status(fig,sprintf('%d points total',size(S.pts,1)));
end

function cbEdit(fig)
% Drag existing points and watch the curvature update live -- the part of the
% CAD workflow that clicking alone cannot reproduce. Needs drawpolyline from
% the Image Processing Toolbox; without it the tool still works, you just
% correct points by undoing and re-clicking instead of nudging them.
S = guidata(fig);
if size(S.pts,1) < 2, status(fig,'add some points first'); return, end
if isempty(which('drawpolyline'))
    status(fig,['dragging needs the Image Processing Toolbox (drawpolyline). ' ...
                'Without it: Undo point and re-click. Everything else works.']);
    return
end
status(fig,'drag any vertex; double-click the line when finished');
roi = drawpolyline(S.axImg,'Position',S.pts,'Color',[1 0.85 0.1], ...
    'LineWidth',1.5,'Closed',S.closed);
l = addlistener(roi,'ROIMoved',@(src,~)onMoved(fig,src));
wait(roi);                       % blocks until the ROI is double-clicked
delete(l);
if isvalid(roi)
    S = guidata(fig); S.pts = roi.Position; guidata(fig,S);
    delete(roi);
end
refit(fig);
end

function onMoved(fig,roi)
% live refit while a vertex is being dragged
S = guidata(fig); S.pts = roi.Position; guidata(fig,S);
S2 = guidata(fig);
if ~isempty(S2.scale) && size(S2.pts,1) >= 4
    refitQuiet(fig);
    F = guidata(fig);
    if ~isempty(F.fit), status(fig,summaryLine(F.fit)); end
end
end

function cbUndo(fig)
S = guidata(fig);
if isempty(S.pts), return, end
S.pts(end,:) = []; guidata(fig,S); refit(fig);
end

function cbClosed(fig)
S = guidata(fig); S.closed = logical(get(S.hClosed,'Value')); guidata(fig,S); refit(fig);
end

function cbSmooth(fig,h)
S = guidata(fig); S.smooth = get(h,'Value'); guidata(fig,S); refit(fig);
end

function cbSave(fig)
S = guidata(fig);
if isempty(S.fit), status(fig,'nothing to save yet'); return, end
[f,p] = uiputfile('*.mat','Save track','mytrack.mat');
if isequal(f,0), return, end
arclength = S.fit.s;
curvature = S.fit.k;
save(fullfile(p,f),'arclength','curvature');
status(fig,sprintf(['saved %s  --  %.1f m, %d stations, %d apexes.\n' ...
    'Put it next to the other tracks in events/ and point Events2 at it.'], ...
    fullfile(p,f),S.fit.L,numel(S.fit.s),S.fit.nApex));
end

function cbLoadTrack(fig)
[f,p] = uigetfile('*.mat','Open an existing track');
if isequal(f,0), return, end
T = load(fullfile(p,f));
if ~all(isfield(T,{'arclength','curvature'}))
    status(fig,'that .mat has no arclength/curvature -- not a track file'); return
end
S = guidata(fig);
S.fit = analyse(T.arclength(:).',T.curvature(:).');
S.pts = zeros(0,2);           % an existing file carries no geometry
guidata(fig,S); redraw(fig);
status(fig,sprintf('inspecting %s (no geometry in the file, curvature only)\n%s', ...
    f,summaryLine(S.fit)));
end

%% ------------------------------------------------------------------- maths
function refit(fig)
computeFit(fig);
redraw(fig);
S = guidata(fig);
if ~isempty(S.fit), status(fig,summaryLine(S.fit)); end
end

function refitQuiet(fig)
% recompute and update only the curvature/radius plots. The image axes are
% left alone because during a drag they hold the live ROI, and redrawing
% them would delete the thing the user is holding.
computeFit(fig);
redrawPlots(fig);
end

function computeFit(fig)
S = guidata(fig);
if size(S.pts,1) < 4 || isempty(S.scale)
    S.fit = []; guidata(fig,S); return
end

% pixels -> metres, y negated so the path sits in a normal right-handed
% frame (image rows increase downward, which would otherwise mirror every
% corner and flip the sign of curvature)
x = S.pts(:,1)*S.scale;
y = -S.pts(:,2)*S.scale;

% Parameterise by cumulative chord length. Duplicated clicks give a zero
% step, which makes the spline singular, so drop them.
keep = [true; hypot(diff(x),diff(y)) > 1e-9];
x = x(keep); y = y(keep);
if numel(x) < 4, S.fit = []; guidata(fig,S); return, end

if S.closed && hypot(x(end)-x(1),y(end)-y(1)) > 1e-9
    x(end+1) = x(1); y(end+1) = y(1);
end
t = [0; cumsum(hypot(diff(x),diff(y)))];

% For a closed loop, wrap a copy of each end onto the other before smoothing
% and keep only the central span afterwards. Smoothing a loop as if it were
% an open line leaves a discontinuity at the seam, which shows up as a corner
% that is not there.
if S.closed
    nw = min(numel(x)-1, max(3,round(0.15*numel(x))));
    xi = [x(end-nw:end-1); x; x(2:nw+1)];
    yi = [y(end-nw:end-1); y; y(2:nw+1)];
    L  = t(end);
    ti = [t(end-nw:end-1)-L; t; t(2:nw+1)+L];
else
    xi = x; yi = y; ti = t;
end

p = smoothParam(mean(diff(ti)),S.smooth);

ppx = csaps(ti,xi,p);
ppy = csaps(ti,yi,p);

% dense resample, then re-parameterise by TRUE arclength: the spline
% parameter is chord length through the clicks, which is not the same thing
tt = linspace(t(1),t(end),20000);
dx = fnval(fnder(ppx,1),tt);  dy = fnval(fnder(ppy,1),tt);
d2x= fnval(fnder(ppx,2),tt);  d2y= fnval(fnder(ppy,2),tt);
sp = hypot(dx,dy);
sArc = [0 cumsum((sp(1:end-1)+sp(2:end))/2 .* diff(tt))];
L = sArc(end);

kk = (dx.*d2y - dy.*d2x)./max(sp.^3,eps);   % signed curvature, 1/m

% output on a uniform arclength grid, first station one step in, matching
% the existing files (Track_Solver prepends the zero itself)
sOut = (1:S.nOut)*L/S.nOut;
kOut = interp1(sArc,kk,sOut,'linear','extrap');

S.fit = analyse(sOut,kOut);
S.fit.xy = [fnval(ppx,tt); fnval(ppy,tt)]/S.scale;   % back to pixels to draw
S.fit.xy(2,:) = -S.fit.xy(2,:);
S.fit.p = p;
guidata(fig,S);
end

function p = smoothParam(h,slider)
% Maps the 0..1 slider onto csaps's smoothing parameter.
%
% p = 1 interpolates every click, including every clicking error; p -> 0
% collapses to a straight line. The useful range sits extremely close to 1
% and scales with the CUBE of the point spacing, so the slider cannot be
% wired to p directly -- the whole usable travel would live in the last 0.1%
% of the bar.
%
% Nor can it be a linear offset from the transition value: 1-p is tiny, so
% subtracting a few multiples of it overshoots straight past zero and the
% slider becomes a cliff between "all the noise" and "a straight line" with
% nothing in between. Measured on the first attempt: 50 apexes at 0.6, one
% apex at 0.4.
%
% Shifting in LOGIT space instead keeps p strictly inside (0,1) however far
% the slider travels, and gives even resolution across the range. Centre is
% MATLAB's own suggested transition value, +/- 3 decades of weight.
pMid = 1/(1+h^3/6);
l = log10(pMid/(1-pMid)) + 6*(slider-0.5);
w = 10^l;
p = w/(1+w);
end

function F = analyse(s,k)
% everything the readout needs, using the SAME apex finder the solver uses
F.s = s(:).'; F.k = k(:).';
F.L = s(end);
F.R = 1./max(abs(k),eps);
try
    e = curvature_apexes(F.s,F.k);
    F.nApex = numel(e);
    r = abs(1./e);
    F.nReal = sum(r < 200);      % apexes tight enough to be corners
    F.minR  = min(r);
catch
    F.nApex = NaN; F.nReal = NaN; F.minR = min(F.R);
end
end

function str = summaryLine(F)
str = sprintf(['length %.1f m   min radius %.2f m   apexes the solver will find: %d ' ...
    '(%d with R < 200 m)\n' ...
    'If those two apex numbers are far apart the curvature is still noisy -- ' ...
    'smooth until they converge on the number of corners you can actually see.'], ...
    F.L,F.minR,F.nApex,F.nReal);
end

%% ------------------------------------------------------------------ drawing
function redraw(fig)
S = guidata(fig);

cla(S.axImg);
if ~isempty(S.img)
    image(S.axImg,S.img); axis(S.axImg,'image'); hold(S.axImg,'on');
end
axis(S.axImg,'off');
if ~isempty(S.pts)
    plot(S.axImg,S.pts(:,1),S.pts(:,2),'o','MarkerSize',4, ...
        'MarkerFaceColor',[1 0.85 0.1],'MarkerEdgeColor','k');
end
if ~isempty(S.fit) && isfield(S.fit,'xy') && ~isempty(S.fit.xy)
    plot(S.axImg,S.fit.xy(1,:),S.fit.xy(2,:),'-','Color',[0.85 0.25 0.1],'LineWidth',2);
end
if isempty(S.scale)
    title(S.axImg,'scale not set');
else
    title(S.axImg,sprintf('%.4f m/px',S.scale));
end

redrawPlots(fig);
end

function redrawPlots(fig)
% curvature and radius only -- safe to call while a drag ROI is live on the
% image axes
S = guidata(fig);
cla(S.axCur); cla(S.axRad);
if isempty(S.fit), return, end

plot(S.axCur,S.fit.s,S.fit.k,'-','LineWidth',1.1);
yline(S.axCur,0,'k:');
xlabel(S.axCur,'distance (m)'); ylabel(S.axCur,'curvature (1/m)');
title(S.axCur,sprintf('%d apexes (%d under 200 m radius)',S.fit.nApex,S.fit.nReal));
grid(S.axCur,'on');

r = min(S.fit.R,200);
plot(S.axRad,S.fit.s,r,'-','LineWidth',1.1);
xlabel(S.axRad,'distance (m)'); ylabel(S.axRad,'radius (m, clipped at 200)');
grid(S.axRad,'on');
end

function status(fig,str)
S = guidata(fig);
set(S.hStat,'String',str);
drawnow limitrate;
end
