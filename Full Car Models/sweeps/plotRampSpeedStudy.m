function fig = plotRampSpeedStudy(results,labels,opts)
%PLOTRAMPSPEEDSTUDY Overlay ramp-speed balance metrics for several cars.
%   fig = plotRampSpeedStudy(results,labels)
%
% RESULTS is a struct array (or cell array) returned by rampSweep, one per
% car. LABELS names those cars. Set opts.outputs to any ordered subset of:
%   aero_front_load, aero_rear_load, aero_balance, mechanical_balance,
%   handling_balance, front_ride_height, rear_ride_height, front_camber,
%   rear_camber, pitch_angle, front_shock_travel, rear_shock_travel.

if nargin < 2 || isempty(labels), labels = strings(0,1); end
if nargin < 3, opts = struct(); end
if iscell(results), results = [results{:}]; end
if isempty(results)
    error('plotRampSpeedStudy:noResults','At least one rampSweep result is required.')
end
if ~isstruct(results)
    error('plotRampSpeedStudy:badResults','results must be rampSweep result structs.')
end

n = numel(results);
labels = string(labels);
labels = labels(:);
if isempty(labels), labels = "car " + string((1:n).'); end
if numel(labels) ~= n
    error('plotRampSpeedStudy:labelCount','labels must contain one name per result.')
end
catalog = outputCatalog();
available = string({catalog.id});
defaultOutputs = ["aero_front_load" "aero_rear_load" "aero_balance" ...
    "mechanical_balance" "handling_balance" "front_ride_height" ...
    "rear_ride_height" "front_camber" "rear_camber" "pitch_angle" ...
    "front_shock_travel" "rear_shock_travel"];
outputs = string(getOr(opts,'outputs',defaultOutputs));
outputs = outputs(:).';
[known,index] = ismember(outputs,available);
if isempty(outputs) || ~all(known)
    bad = strjoin(outputs(~known),', ');
    error('plotRampSpeedStudy:unknownOutput', ...
        'Unknown output(s): %s. Choose from: %s.',bad,strjoin(available,', '))
end
selected = catalog(index);
for i = 1:n
    if ~isfield(results(i),'perSpeed') || ~istable(results(i).perSpeed)
        error('plotRampSpeedStudy:badResult','Result %d has no perSpeed table.',i)
    end
    required = [{'vCar'} {selected.field}];
    if ~all(ismember(required,results(i).perSpeed.Properties.VariableNames))
        error('plotRampSpeedStudy:missingMetrics', ...
            'Result %d is missing a metric required by the selected outputs.',i)
    end
end

visible = getOr(opts,'visible','on');
name = getOr(opts,'figureName','ramp speed study: balance comparison');
columns = min(3,numel(selected));
rows = ceil(numel(selected)/columns);
fig = figure('Name',name,'Visible',visible, ...
    'Position',[60 60 430*columns 350*rows]);
setPlotFont(fig);
tiledlayout(fig,rows,columns,'TileSpacing','compact','Padding','compact');
cmap = lines(max(n,2));

axesHandles = gobjects(numel(selected),1);
for j = 1:numel(selected)
    spec = selected(j);
    ax = nexttile; axesHandles(j) = ax; hold(ax,'on');
    if spec.zeroLine
        yline(ax,0,'k:','HandleVisibility','off','LineWidth',1.1);
    end
    for i = 1:n
        S = results(i).perSpeed;
        plot(ax,S.vCar,spec.scale*S.(spec.field),'-o','Color',cmap(i,:), ...
            'LineWidth',1.8,'MarkerSize',4,'DisplayName',labels(i));
    end
    formatAxes(ax,'vCar (m/s)',spec.yLabel,spec.title);
    if strlength(spec.subtitle) > 0
        subtitle(ax,spec.subtitle,'FontSize',8);
    end
    legend(ax,'Location','best','FontSize',7);
end
linkaxes(axesHandles,'x');
sgtitle(fig,'ramp-speed comparison: aero, ride height, camber, and handling balance');
end

function formatAxes(ax,xText,yText,titleText)
grid(ax,'on'); box(ax,'on');
xlabel(ax,xText); ylabel(ax,yText); title(ax,titleText);
end

function catalog = outputCatalog()
catalog = struct( ...
    'id',{'aero_front_load','aero_rear_load','aero_balance','mechanical_balance', ...
          'handling_balance','front_ride_height','rear_ride_height', ...
          'front_camber','rear_camber','pitch_angle', ...
          'front_shock_travel','rear_shock_travel'}, ...
    'field',{'aero_downforce_front_N','aero_downforce_rear_N','aero_balance','mech_balance', ...
             'K_linear','front_ride_height_in','rear_ride_height_in', ...
             'front_camber_deg','rear_camber_deg','pitch_angle_deg', ...
             'front_shock_travel_in','rear_shock_travel_in'}, ...
    'scale',{1,1,100,100,1,1,1,1,1,1,1,1}, ...
    'yLabel',{'front aero load (N)','rear aero load (N)','front aero balance (%)', ...
              'front mechanical balance (%)','K_{linear} (deg/g)', ...
              'front ride height (in)','rear ride height (in)', ...
              'front camber magnitude (deg)','rear camber magnitude (deg)', ...
              'pitch angle (deg)','front shock compression (in)', ...
              'rear shock compression (in)'}, ...
    'title',{'front aerodynamic load','rear aerodynamic load','centre of pressure', ...
             'mechanical balance (LLTD)','total handling balance', ...
             'front ride height','rear ride height','front camber','rear camber', ...
             'vehicle pitch','front shock travel','rear shock travel'}, ...
    'subtitle',{'','','','','positive = understeer, negative = oversteer','','','','', ...
                'positive = nose-up, negative = nose-down','positive = compression', ...
                'positive = compression'}, ...
    'zeroLine',{false,false,false,false,true,false,false,false,false,true,false,false});
end

function value = getOr(s,name,defaultValue)
if isfield(s,name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end
end
