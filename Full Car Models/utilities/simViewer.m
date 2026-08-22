function simViewer(src,event,car)
% Channel viewer for solved lap data -- add and remove channels, scrub a
% cursor, and see the track map coloured by whatever channel is selected.
%
%   simViewer(ev)                 an Events2 object that has been run
%   simViewer(ev,'endurance')     pick the event (default 'autocross')
%   simViewer(ev,'autocross',car) also derive envelope-based channels
%   simViewer(traceStruct)        a struct with t/v/gLat/gLong (from aeroLapTimes)
%   simViewer                     file picker for a .mat holding either
%
%
% DERIVED CHANNELS ARE EXACT, NOT MODELLED
% yaw rate is gLat*g/v, path radius is v^2/(gLat*g), theta is the g-g polar
% angle. These are identities, not extra physics. The one channel that pulls
% in outside information is "grip used", which needs the car's own g-g -- it
% is the fraction of the available envelope the lap is actually using at each
% point, and it only appears when a solved car is supplied.
%


if nargin < 1 || isempty(src), src = pickFile(); if isempty(src), return, end, end
if nargin < 2 || isempty(event), event = 'autocross'; end
if nargin < 3, car = []; end

C = buildChannels(src,event,car);
if isempty(C.chan)
    error('simViewer:noData','no usable channels -- has the event been run?');
end

S = struct();
S.C        = C;
S.sel      = defaultSelection(C);
S.xmode    = 1;                      % 1 = distance, 2 = time
S.cursor   = C.x1(1);
S.colorBy  = firstSel(S.sel);
% initialised here so every handler can test them before anything is drawn
S.map      = [];
S.hDot     = gobjects(1);
S.hCur     = gobjects(0);
S.ax       = gobjects(0);

S.fig = figure('Name',sprintf('simViewer - %s',C.event),'NumberTitle','off', ...
    'Position',[60 60 1400 800],'Color',[0.94 0.94 0.94]);
setPlotFont(S.fig);

uicontrol(S.fig,'Style','text','String','channels', ...
    'Units','normalized','Position',[0.012 0.955 0.13 0.03], ...
    'BackgroundColor',[0.94 0.94 0.94],'HorizontalAlignment','left','FontSize',9);
S.hList = uicontrol(S.fig,'Style','listbox','Max',10,'Min',0, ...
    'String',{C.chan.label},'Value',find(S.sel), ...
    'Units','normalized','Position',[0.012 0.34 0.15 0.615], ...
    'Callback',@(h,~)cbSelect(S.fig,h));

uicontrol(S.fig,'Style','text','String', ...
    'ctrl-click to add/remove.  * = solved, no mark = exact algebra', ...
    'Units','normalized','Position',[0.012 0.295 0.15 0.04], ...
    'BackgroundColor',[0.94 0.94 0.94],'HorizontalAlignment','left','FontSize',7);

S.hX = uicontrol(S.fig,'Style','popupmenu','String',{'x: distance (m)','x: time (s)'}, ...
    'Units','normalized','Position',[0.012 0.245 0.15 0.04], ...
    'Callback',@(h,~)cbX(S.fig,h));

uicontrol(S.fig,'Style','pushbutton','String','Export CSV', ...
    'Units','normalized','Position',[0.012 0.20 0.07 0.038], ...
    'Callback',@(~,~)cbExport(S.fig));
uicontrol(S.fig,'Style','pushbutton','String','Select all', ...
    'Units','normalized','Position',[0.092 0.20 0.07 0.038], ...
    'Callback',@(~,~)cbAll(S.fig));

S.hRead = uicontrol(S.fig,'Style','text','String','', ...
    'Units','normalized','Position',[0.012 0.012 0.15 0.18], ...
    'BackgroundColor',[1 1 1],'HorizontalAlignment','left','FontSize',8);

S.axMap = axes('Parent',S.fig,'Units','normalized','Position',[0.185 0.06 0.24 0.34]);
S.panel = uipanel('Parent',S.fig,'Units','normalized','Position',[0.185 0.50 0.80 0.48], ...
    'BackgroundColor',[1 1 1],'BorderType','none');

% Scrub bar. Spans the same range as the trace x-axis, so dragging it walks
% the dot around the track map. Clicking directly on a trace also works and
% moves the slider to match -- they are two views of one cursor.
uicontrol(S.fig,'Style','text','String','lap position', ...
    'Units','normalized','Position',[0.185 0.435 0.08 0.03], ...
    'BackgroundColor',[0.94 0.94 0.94],'HorizontalAlignment','left','FontSize',8);
S.hSlider = uicontrol(S.fig,'Style','slider', ...
    'Min',C.x1(1),'Max',C.x1(end),'Value',C.x1(1), ...
    'Units','normalized','Position',[0.185 0.415 0.80 0.025], ...
    'Callback',@(h,~)cbSlider(S.fig,h));
% continuous drag rather than only on release
try
    addlistener(S.hSlider,'ContinuousValueChange',@(h,~)cbSlider(S.fig,h));
catch
    % older releases without the listener: the slider still works on release
end

uicontrol(S.fig,'Style','text','String','colour map by', ...
    'Units','normalized','Position',[0.46 0.435 0.08 0.03], ...
    'BackgroundColor',[0.94 0.94 0.94],'HorizontalAlignment','left','FontSize',8);
S.hColor = uicontrol(S.fig,'Style','popupmenu','String',{C.chan.label}, ...
    'Value',S.colorBy,'Units','normalized','Position',[0.54 0.432 0.16 0.035], ...
    'Callback',@(h,~)cbColor(S.fig,h));

guidata(S.fig,S);
rebuild(S.fig);
end

%% ------------------------------------------------------------- channel build
function C = buildChannels(src,event,car)
C.event = event;
C.chan  = struct('label',{},'name',{},'unit',{},'data',{},'solved',{});

[t,v,ay,ax,s,kappa] = extract(src,event);
if isempty(v), C.x1 = []; C.x2 = []; return, end

n = numel(v);
if isempty(s)
    % no track paired with this trace -- integrate velocity for distance so
    % the distance axis still means something
    s = [0 cumsum((v(1:end-1)+v(2:end))/2 .* diff(t(:).'))];
    kappa = nan(1,n);
end

g = 9.81;
C.x1 = s(:).';  C.x2 = t(:).';        % distance and time axes

add = @(lab,nm,u,d,sv) struct('label',lab,'name',nm,'unit',u,'data',d(:).','solved',sv);

C.chan(end+1) = add('* vCar (m/s)','vCar','m/s',v,true);
C.chan(end+1) = add('  vCar (km/h)','vCar_kph','km/h',v*3.6,false);
C.chan(end+1) = add('* gLat (g)','gLat','g',ay/g,true);
C.chan(end+1) = add('* gLong (g)','gLong','g',ax/g,true);

gLat = ay(:).'/g; gLong = ax(:).'/g;
C.chan(end+1) = add('  gMag (g)','gMag','g',hypot(gLat,gLong),false);
C.chan(end+1) = add('  theta (deg)','theta','deg',atan2d(abs(gLat),gLong),false);

% identities, not extra physics
yaw = gLat*g./max(v(:).',eps);
C.chan(end+1) = add('  yaw rate (deg/s)','yaw_rate','deg/s',yaw*180/pi,false);
R = v(:).'.^2./max(abs(gLat)*g,eps);
C.chan(end+1) = add('  path radius (m)','path_radius','m',min(R,500),false);

if any(isfinite(kappa))
    C.chan(end+1) = add('* track curvature (1/m)','curvature','1/m',kappa,true);
    C.chan(end+1) = add('  track radius (m)','track_radius','m', ...
        min(1./max(abs(kappa),eps),500),false);
end

% longitudinal power crossing the contact patches: m*a*v. Not engine output
% -- it excludes drag and drivetrain loss -- but it is exactly what the solved
% acceleration implies.
if ~isempty(car) && isprop(car,'M')
    C.chan(end+1) = add('  long. power (kW)','power','kW',car.M*ax(:).'.*v(:).'/1000,false);
end

% envelope usage: how much of the g-g the lap is actually leaning on. The one
% channel that needs information from outside the trace.
if ~isempty(car) && isprop(car,'ss_info') && ~isempty(car.ss_info)
    u = gripUsage(car,v(:).',gLat,gLong);
    if any(isfinite(u))
        C.chan(end+1) = add('  grip used (-)','grip_used','-',u,false);
    end
end
end

function [t,v,ay,ax,s,kappa] = extract(src,event)
t=[]; v=[]; ay=[]; ax=[]; s=[]; kappa=[]; %#ok<NASGU>
if isa(src,'Events2')
    if ~isprop(src,event) || isempty(src.(event))
        error('simViewer:notRun','%s has not been run on that Events2 object',event);
    end
    E = src.(event);
    if ~isfield(E,'long_vel'), error('simViewer:noTrace','%s holds no trace',event); end
    v = E.long_vel(:).'; t = E.time_vec(:).';
    if isfield(E,'lat_accel'),  ay = E.lat_accel(:).';  end
    if isfield(E,'long_accel'), ax = E.long_accel(:).'; end
    trk = [];
    if strcmpi(event,'autocross'), trk = src.autocross_track;
    elseif strcmpi(event,'endurance'), trk = src.endurance_track; end
    if ~isempty(trk), s = trk(1,:); kappa = trk(2,:); end
elseif isstruct(src) && isfield(src,'v')
    v = src.v(:).'; t = src.t(:).';
    if isfield(src,'gLat'),  ay = src.gLat(:).'*9.81;  end
    if isfield(src,'gLong'), ax = src.gLong(:).'*9.81; end
else
    error('simViewer:badInput','give an Events2 object or a trace struct with t/v');
end

% these vectors routinely differ by a sample or two depending on how the
% solver terminated; trim rather than guess which end is stale
m = min([numel(t) numel(v) numel(ay) numel(ax)]);
if ~isempty(s), m = min(m,numel(s)); s = s(1:m); kappa = kappa(1:m); end
t=t(1:m); v=v(1:m); ay=ay(1:m); ax=ax(1:m);
end

function u = gripUsage(car,v,gLat,gLong)
% fraction of the solved g-g envelope in use at each point
u = nan(size(v));
try
    grids.vCar  = unique(round(car.ss_info(:,6),6));
    grids.theta = linspace(0,180,73);
    grids.gLat  = linspace(0,2.5,41);
    grids.gLong = linspace(-3,2,41);
    E = aeroEnvelope(car,grids);
    th = atan2d(abs(gLat),gLong);
    r0 = interp2(grids.theta,grids.vCar(:),E.r_theta,th,v,'linear',NaN);
    u  = hypot(gLat,gLong)./r0;
catch
    % no envelope, no channel -- caller drops it
end
end

%%  GUI
function rebuild(fig)
S = guidata(fig);
delete(allchild(S.panel));
idx = find(S.sel);
if isempty(idx)
    uicontrol(S.panel,'Style','text','String','no channels selected', ...
        'Units','normalized','Position',[0.4 0.48 0.2 0.06],'BackgroundColor',[1 1 1]);
    drawMap(fig); return
end

x = xdata(S);
n = numel(idx);
S.ax = gobjects(1,n);
h = 0.98/n;
for k = 1:n
    c = S.C.chan(idx(k));
    ax = axes('Parent',S.panel,'Units','normalized', ...
        'Position',[0.075 0.99-k*h 0.90 h*0.86]);
    plot(ax,x,c.data,'-','LineWidth',1.1,'Color',[0.10 0.35 0.70]);
    grid(ax,'on'); box(ax,'on');
    ylabel(ax,sprintf('%s (%s)',c.name,c.unit),'FontSize',8,'Interpreter','none');
    set(ax,'FontSize',8);
    xlim(ax,[min(x) max(x)]);
    if k < n, set(ax,'XTickLabel',[]); end
    S.ax(k) = ax;
end
xlabel(S.ax(end),xlabelStr(S),'FontSize',9);
% one cursor line per axes, all moved together
S.hCur = gobjects(1,n);
for k = 1:n
    S.hCur(k) = xline(S.ax(k),S.cursor,'-','Color',[0.85 0.25 0.1],'LineWidth',1.2);
end
set(S.fig,'WindowButtonDownFcn',@(f,~)onClick(f));
guidata(fig,S);
drawMap(fig);
updateReadout(fig);
end

function drawMap(fig)


S = guidata(fig);
cla(S.axMap); S.hDot = gobjects(1);
k = findChan(S.C,'curvature');
if isempty(k)
    text(S.axMap,0.5,0.5,'no track curvature, so no map', ...
        'HorizontalAlignment','center','Parent',S.axMap);
    axis(S.axMap,'off'); S.map = []; guidata(fig,S); return
end
s = S.C.x1; kap = S.C.chan(k).data;
ok = isfinite(s) & isfinite(kap);
sOK = s(ok); kOK = kap(ok);
psi = cumtrapz(sOK,kOK);
X = cumtrapz(sOK,cos(psi));  Y = cumtrapz(sOK,sin(psi));

% decimate for drawing only; the dot still indexes the full-resolution path
step = max(1,round(numel(X)/4000));
col  = S.C.chan(S.colorBy).data(ok);
scatter(S.axMap,X(1:step:end),Y(1:step:end),6,col(1:step:end),'filled');
axis(S.axMap,'equal'); axis(S.axMap,'off');
cb = colorbar(S.axMap);
cb.Label.String = S.C.chan(S.colorBy).name; cb.Label.Interpreter = 'none';
title(S.axMap,sprintf('track coloured by %s',S.C.chan(S.colorBy).name), ...
    'Interpreter','none','FontSize',9);
hold(S.axMap,'on');
S.hDot = plot(S.axMap,X(1),Y(1),'o','MarkerSize',10,'LineWidth',2, ...
    'MarkerEdgeColor',[0.85 0.25 0.1],'MarkerFaceColor',[1 0.95 0.9]);

S.map = struct('X',X,'Y',Y,'sMap',sOK,'tMap',S.C.x2(ok));
guidata(fig,S);
moveDot(fig);
end

function moveDot(fig)
% cheap: index the prebuilt path and move one marker
S = guidata(fig);
if isempty(S.map) || ~isgraphics(S.hDot), return, end
if S.xmode == 1, ref = S.map.sMap; else, ref = S.map.tMap; end
[~,i] = min(abs(ref - S.cursor));
set(S.hDot,'XData',S.map.X(i),'YData',S.map.Y(i));
end

function moveCursor(fig,xval)
% one place that moves the cursor, whatever drove it: slider or a click on a
% trace. Keeps the slider, the trace lines, the map dot and the readout from
% ever disagreeing about where "now" is.
S = guidata(fig);
x = xdata(S);
S.cursor = min(max(xval,min(x)),max(x));
guidata(fig,S);
for k = 1:numel(S.hCur)
    if isgraphics(S.hCur(k)), S.hCur(k).Value = S.cursor; end
end
if isgraphics(S.hSlider), set(S.hSlider,'Value',S.cursor); end
moveDot(fig);
updateReadout(fig);
end

function onClick(fig)
S = guidata(fig);
if ~isfield(S,'ax') || isempty(S.ax), return, end
p = get(S.ax(1),'CurrentPoint');
inAny = false;
for k = 1:numel(S.ax)
    pk = get(S.ax(k),'CurrentPoint');
    yl = get(S.ax(k),'YLim'); xl = get(S.ax(k),'XLim');
    if pk(1,1)>=xl(1) && pk(1,1)<=xl(2) && pk(1,2)>=yl(1) && pk(1,2)<=yl(2)
        p = pk; inAny = true; break
    end
end
if ~inAny, return, end
moveCursor(fig,p(1,1));
end

function cbSlider(fig,h)
moveCursor(fig,get(h,'Value'));
end

function cbColor(fig,h)
S = guidata(fig); S.colorBy = get(h,'Value'); guidata(fig,S);
drawMap(fig);
end

function updateReadout(fig)
S = guidata(fig);
x = xdata(S);
[~,i] = min(abs(x-S.cursor));
str = sprintf('%s = %.2f\n\n',xlabelStr(S),x(i));
for k = find(S.sel)
    c = S.C.chan(k);
    str = [str sprintf('%-14s %9.3f %s\n',c.name,c.data(i),c.unit)]; %#ok<AGROW>
end
set(S.hRead,'String',str);
end

%% -------------------------------------------------------------- callbacks
function cbSelect(fig,h)
S = guidata(fig);
S.sel(:) = false; S.sel(get(h,'Value')) = true;
if ~S.sel(S.colorBy), S.colorBy = firstSel(S.sel); end
guidata(fig,S); rebuild(fig);
end

function cbX(fig,h)
S = guidata(fig); S.xmode = get(h,'Value');
% the slider spans whichever axis is showing, so its limits move with it --
% leaving them on the old range would silently clamp the cursor
x = xdata(S);
S.cursor = x(1);
set(S.hSlider,'Min',min(x),'Max',max(x),'Value',x(1));
guidata(fig,S); rebuild(fig);
end

function cbAll(fig)
S = guidata(fig); S.sel(:) = true;
set(S.hList,'Value',1:numel(S.sel));
guidata(fig,S); rebuild(fig);
end

function cbExport(fig)
S = guidata(fig);
[f,p] = uiputfile('*.csv','Export channels','lap_channels.csv');
if isequal(f,0), return, end
T = table(S.C.x1(:),S.C.x2(:),'VariableNames',{'distance_m','time_s'});
for k = find(S.sel)
    T.(matlab.lang.makeValidName(S.C.chan(k).name)) = S.C.chan(k).data(:);
end
writetable(T,fullfile(p,f));
set(S.hRead,'String',sprintf('wrote %s\n%d rows, %d channels', ...
    fullfile(p,f),height(T),width(T)));
end

%% ----------------------------------------------------------------- helpers
function x = xdata(S)
if S.xmode == 1, x = S.C.x1; else, x = S.C.x2; end
end
function s = xlabelStr(S)
if S.xmode == 1, s = 'distance (m)'; else, s = 'time (s)'; end
end
function k = findChan(C,name)
k = find(strcmp({C.chan.name},name),1);
end
function i = firstSel(sel)
i = find(sel,1); if isempty(i), i = 1; end
end
function sel = defaultSelection(C)
sel = false(1,numel(C.chan));
for nm = {'vCar','gLat','gLong'}
    k = find(strcmp({C.chan.name},nm{1}),1);
    if ~isempty(k), sel(k) = true; end
end
if ~any(sel), sel(1:min(3,numel(sel))) = true; end
end
function src = pickFile()
src = [];
[f,p] = uigetfile('*.mat','Open solved lap data');
if isequal(f,0), return, end
D = load(fullfile(p,f));
fn = fieldnames(D);
for i = 1:numel(fn)
    v = D.(fn{i});
    if isa(v,'Events2'), src = v; return, end
    if iscell(v) && ~isempty(v) && isstruct(v{1}) && isfield(v{1},'autocross')
        src = v{1}.autocross; return          % traces cell from aeroLapTimes
    end
    if isstruct(v) && isfield(v,'v') && isfield(v,'t'), src = v; return, end
end
error('simViewer:nothingUsable','%s holds no Events2 object or trace struct',f);
end
