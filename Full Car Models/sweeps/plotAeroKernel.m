function plotAeroKernel(S,event,param)
% Time-at-condition kernels from the solved lap traces: where the car
% actually spends its lap. Read alongside the sensitivity plots these say
% whether a gain sits somewhere that matters.
%
% inputs  S:     output of aeroSensitivity (must have been given traces)
%         event: 'autocross' (default) | 'endurance'
%         param: 'ClA' | 'CdA' | 'CoP' -- draw the baseline against just this
%                parameter's fastest case. Omit to draw the baseline against
%                the fastest of every swept parameter at once.
%
% The per-parameter form is the one to use when sweeping all three: three
% fastest cases on one axis obscures the very differences the figure exists
% to show, since the kernels sit almost on top of each other.

if nargin < 2 || isempty(event), event = 'autocross'; end
if nargin < 3, param = ''; end
if ~isfield(S,'traces') || isempty(S.traces)
    error('plotAeroKernel:noTraces', ...
        ['no lap traces in S -- run aeroLapTimes with two outputs and pass ' ...
         'them in as opts.traces']);
end

base = S.baseIdx;

% which cases to draw: baseline plus the fastest of the requested parameter,
% or of every swept parameter when none was named
if isempty(param)
    pf = fieldnames(S.sens);
else
    if ~isfield(S.sens,param)
        error('plotAeroKernel:noParam','no sensitivity computed for %s',param);
    end
    pf = {param};
end

show = base;
names = {'baseline'};
for i = 1:numel(pf)
    if isfield(S.sens.(pf{i}),'fastestCase')
        fc = S.sens.(pf{i}).fastestCase;
        if fc ~= base && ~ismember(fc,show)
            show(end+1) = fc; %#ok<AGROW>
            names{end+1} = sprintf('fastest %s (case %d)',pf{i},fc); %#ok<AGROW>
        end
    end
end

vars = {'vCar','v (m/s)'; 'gLat','|gLat| (g)'; 'gLong','gLong (g)'; 'theta','\theta (deg)'};

tag = '';
if ~isempty(param), tag = sprintf(' - %s',param); end
figure('Name',sprintf('Lap kernels (%s)%s',event,tag),'Position',[80 80 1150 720]);
setPlotFont(gcf);   % must precede the subplots: these are DEFAULT properties
co = lines(numel(show));
for k = 1:size(vars,1)
    subplot(2,2,k); hold on
    for j = 1:numel(show)
        tr = getTrace(S.traces{show(j)},event);
        [d,g,info] = lapKernel(tr,vars{k,1});
        if all(isnan(d)), continue, end
        plot(g,d,'-','Color',co(j,:),'LineWidth',2,'DisplayName',names{j});
        if j == 1
            % shade the baseline so the eye reads it as the weighting
            fill([g(:);flipud(g(:))],[d(:);zeros(numel(d),1)], co(j,:), ...
                'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
            xline(info.mean,'--','Color',co(j,:),'HandleVisibility','off');
        end
    end
    grid on; box on;
    xlabel(vars{k,2}); ylabel('time density (1/unit)');
    title(sprintf('time spent vs %s',vars{k,1}));
    if k == 1, legend('Location','best','FontSize',8); end
end
sgtitle(sprintf('Where the lap is actually spent -- %s%s (dashed = time-weighted mean)', ...
    event, tag));
end

function tr = getTrace(t,event)
if isempty(t) || ~isfield(t,event) || isempty(t.(event))
    tr = [];
else
    tr = t.(event);
end
end
