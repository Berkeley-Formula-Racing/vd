function tests = test_rampSpeedPlotCatalog
tests = functiontests(localfunctions);
end

function testCatalogCoversPlotMetricFamilies(testCase)
catalog = rampSpeed.metricCatalog();
ids = string({catalog.id});
expected = ["aero_front_load","aero_rear_load","aero_balance", ...
    "mechanical_balance","handling_balance","yaw_balance", ...
    "capability_free","capability_sustainable","ramp_complete", ...
    "downforce","drag","front_ride_height","rear_ride_height", ...
    "front_camber","rear_camber","pitch_angle", ...
    "front_shock_travel","rear_shock_travel","front_axle_load", ...
    "rear_axle_load","min_wheel_load","wheel_lift","validity", ...
    "truncated","power_limited","raw_aLat","raw_aLong","raw_steer"];
verifyTrue(testCase,all(ismember(expected,ids)));
requiredFields = ["id","tab","validTypes","sourceLevel","field", ...
    "derivation","storedUnits","displayUnits","scale","yLabel", ...
    "title","subtitle","zeroLine","seriesStyle","validityRule"];
verifyTrue(testCase,all(isfield(catalog,cellstr(requiredFields))));
longitudinal = rampSpeed.metricCatalog(struct("type","longitudinal"));
verifyFalse(testCase,ismember("capability_free",string({longitudinal.id})));
verifyTrue(testCase,ismember("drag",string({longitudinal.id})));
end

function testCatalogPreservesLegacyMetricContract(testCase)
catalog = rampSpeed.metricCatalog();
spec = catalog(strcmp(string({catalog.id}),"aero_front_load"));
verifyEqual(testCase,spec.sourceLevel,"perSpeed");
verifyEqual(testCase,spec.field,"front_downforce_N");
verifyEqual(testCase,spec.storedUnits,"N");
verifyEqual(testCase,spec.displayUnits,"N");
verifyEqual(testCase,spec.scale,1);
verifyFalse(testCase,spec.zeroLine);
handling = catalog(strcmp(string({catalog.id}),"handling_balance"));
verifyTrue(testCase,isa(handling.derivation,"function_handle"));
verifyEqual(testCase,handling.displayUnits,"deg/g");
end

function testBuildPlotDataUsesNativeGridAndCarriesMasks(testCase)
runs(1) = makePlotRun("baseline",[5;10;15],[100;200;300], ...
    [true;false;true],[false;false;true]);
runs(2) = makePlotRun("variant",[5;12.5;15],[90;190;290], ...
    [true;true;true],[false;false;false]);
data = rampSpeed.buildPlotData(runs,"aero_front_load", ...
    struct("labels",["base";"variant"]));
verifyEqual(testCase,numel(data.series),2);
verifyEqual(testCase,data.series(1).x,[5;10;15]);
verifyEqual(testCase,data.series(2).x,[5;12.5;15]);
verifyEqual(testCase,data.series(1).values(1),100);
verifyTrue(testCase,isnan(data.series(1).values(2)));
verifyEqual(testCase,data.series(1).valid,[true;false;true]);
verifyEqual(testCase,data.series(1).truncated,[false;false;true]);
verifyEqual(testCase,data.series(1).power_limited,[false;false;false]);
verifyEqual(testCase,data.series(1).wheel_lift,[false;false;false]);
verifyEqual(testCase,data.series(1).reason(2),"invalid speed");
end

function testBuildPlotDataGroupsRawPointsBySpeedIndex(testCase)
run = makePlotRun("baseline",[5;10],[100;200],[true;true],[false;false]);
run.points = makePointTable(run.points,4);
run.points.speed_mps = [5;5;5;5];
run.points.speed_index = [1;1;2;2];
run.points.point_index = [1;2;1;2];
run.points.aLat_mps2 = [1;1.1;2;2.1]*9.80665;
run.points.valid = true(4,1);
run.points.status(:) = "complete";
data = rampSpeed.buildPlotData(run,"raw_aLat");
verifyEqual(testCase,data.sourceLevel,"points");
verifyEqual(testCase,numel(data.series(1).groups),2);
verifyEqual(testCase,[data.series(1).groups.speedIndex],[1 2]);
verifyEqual(testCase,data.series(1).groups(1).x,[5;5]);
verifyEqual(testCase,data.series(1).groups(2).values,[2;2.1]);
end

function testBuildInspectorTableIncludesSIAndDisplayColumns(testCase)
run = makePlotRun("baseline",5,100,true,false);
selection = struct("speed_mps",5,"speedIndex",1,"pointIndex",1);
units = struct("speed","mph","force","lbf","length","in", ...
    "angle","deg","acceleration","g","angularRate","rps");
T = rampSpeed.buildInspectorTable(run,selection,units);
required = ["speed_mps","speed_display","Fz_FL_N","Fz_FL_display", ...
    "Fx_FR_N","Fx_FR_display","gamma_RL_rad","gamma_RL_display", ...
    "valid","truncated","power_limited","wheel_lift", ...
    "aero_outside_map","max_constraint_residual"];
verifyTrue(testCase,all(ismember(required,string(T.Properties.VariableNames))));
verifyEqual(testCase,height(T),1);
verifyEqual(testCase,T.speed_display,5*2.2369362920544,"AbsTol",1e-12);
verifyEqual(testCase,T.Fz_FL_display,T.Fz_FL_N/4.4482216152605, ...
    "AbsTol",1e-12);
verifyEqual(testCase,T.gamma_RL_display,T.gamma_RL_rad*180/pi, ...
    "AbsTol",1e-12);
end

function testRenderMetricDrawsWarningsAndZeroLine(testCase)
run = makePlotRun("baseline",[5;10;15],[100;200;300], ...
    [true;true;true],[false;true;false]);
data = rampSpeed.buildPlotData(run,"handling_balance");
fig = figure("Visible","off");
cleanup = onCleanup(@()close(fig));
ax = axes(fig);
h = rampSpeed.renderMetric(ax,data,struct("showLegend",false));
verifyTrue(testCase,isfield(h,"lines"));
verifyGreaterThanOrEqual(testCase,numel(h.lines),1);
verifyTrue(testCase,isfield(h,"warningMarkers"));
verifyGreaterThanOrEqual(testCase,numel(h.warningMarkers),1);
verifyEqual(testCase,numel(findall(ax,"Type","constantline")),1);
verifyEqual(testCase,string(ax.YLabel.String),data.yLabel);
end

function run = makePlotRun(id,speeds,frontLoads,valid,truncated)
speeds = double(speeds(:));
n = numel(speeds);
if isscalar(frontLoads)
    frontLoads = repmat(double(frontLoads),n,1);
else
    frontLoads = double(frontLoads(:));
end
valid = logical(valid);
valid = valid(:);
truncated = logical(truncated);
truncated = truncated(:);
caseInfo = struct("id",string(id),"label",string(id),"source","test", ...
    "designRow",1,"sourceIndex",1,"carRole","lap","carColumn",1);
settings = struct("speeds",speeds.');
run = rampSpeed.makeRun("lateral","coast",settings,caseInfo);
run.perSpeed.front_downforce_N = frontLoads;
run.perSpeed.rear_downforce_N = 0.75*frontLoads;
run.perSpeed.downforce_N = 1.75*frontLoads;
run.perSpeed.drag_N = 0.1*frontLoads;
run.perSpeed.aero_balance_front = 0.48 + 0.01*(1:n).';
run.perSpeed.mechanical_balance_front = 0.52 + 0.01*(1:n).';
run.perSpeed.K_linear_rad_per_mps2 = (1:n).'/(180/pi*9.80665);
run.perSpeed.cuo_steer_linear_rad = (1:n).'/1000;
run.perSpeed.aLat_free_mps2 = 1.5*9.80665*ones(n,1);
run.perSpeed.aLat_sustainable_mps2 = 1.2*9.80665*ones(n,1);
run.perSpeed.ramp_complete_fraction = 1 - 0.1*double(truncated);
run.perSpeed.front_ride_height_m = 0.03*ones(n,1);
run.perSpeed.rear_ride_height_m = 0.025*ones(n,1);
run.perSpeed.front_camber_rad = 0.01*ones(n,1);
run.perSpeed.rear_camber_rad = -0.012*ones(n,1);
run.perSpeed.pitch_rad = 0.002*ones(n,1);
run.perSpeed.front_shock_travel_m = 0.002*ones(n,1);
run.perSpeed.rear_shock_travel_m = 0.001*ones(n,1);
run.perSpeed.Fz_front_axle_N = 900*ones(n,1);
run.perSpeed.Fz_rear_axle_N = 700*ones(n,1);
run.perSpeed.min_Fz_N = 150*ones(n,1);
run.perSpeed.valid = valid;
run.perSpeed.status(:) = "complete";
run.perSpeed.reason(:) = "";
run.perSpeed.reason(~valid) = "invalid speed";
run.perSpeed.truncated = truncated;
run.perSpeed.power_limited = false(n,1);
run.perSpeed.wheel_lift = false(n,1);
run.perSpeed.aero_outside_map = false(n,1);
run.perSpeed.aero_residual_m = 0.001*ones(n,1);
run.status = "complete";
run.runMeta.status = "complete";
run.points = makePointTable(run.points,n);
run.points.speed_mps = speeds;
run.points.speed_index = (1:n).';
run.points.point_index = ones(n,1);
run.points.valid = valid;
run.points.status(:) = "complete";
run.points.aLat_mps2 = 1.1*9.80665*ones(n,1);
run.points.aLong_mps2 = zeros(n,1);
run.points.Fz_FL_N = 450*ones(n,1);
run.points.Fz_FR_N = 440*ones(n,1);
run.points.Fz_RL_N = 430*ones(n,1);
run.points.Fz_RR_N = 425*ones(n,1);
run.points.Fx_FL_N = 100*ones(n,1);
run.points.Fx_FR_N = 101*ones(n,1);
run.points.Fx_RL_N = 90*ones(n,1);
run.points.Fx_RR_N = 91*ones(n,1);
run.points.Fy_FL_N = 220*ones(n,1);
run.points.Fy_FR_N = 221*ones(n,1);
run.points.Fy_RL_N = 210*ones(n,1);
run.points.Fy_RR_N = 211*ones(n,1);
run.points.gamma_FL_rad = 0.01*ones(n,1);
run.points.gamma_FR_rad = 0.011*ones(n,1);
run.points.gamma_RL_rad = 0.012*ones(n,1);
run.points.gamma_RR_rad = 0.013*ones(n,1);
run.points.alpha_FL_rad = 0.02*ones(n,1);
run.points.alpha_FR_rad = 0.021*ones(n,1);
run.points.alpha_RL_rad = 0.022*ones(n,1);
run.points.alpha_RR_rad = 0.023*ones(n,1);
run.points.kappa_FL = 0.01*ones(n,1);
run.points.kappa_FR = 0.011*ones(n,1);
run.points.kappa_RL = 0.012*ones(n,1);
run.points.kappa_RR = 0.013*ones(n,1);
run.points.aero_outside_map = false(n,1);
run.points.aero_residual_m = 0.001*ones(n,1);
run.points.max_constraint_residual = 0.01*ones(n,1);
run.points.max_equality_residual = 0.005*ones(n,1);
run.points.max_inequality_violation = 0.002*ones(n,1);
end

function T = makePointTable(template,n)
variableNames = template.Properties.VariableNames;
variableTypes = varfun(@class,template,'OutputFormat','cell');
T = table('Size',[n width(template)], ...
    'VariableTypes',variableTypes,'VariableNames',variableNames);
end
