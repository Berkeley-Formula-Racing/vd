function [cars,cases,designTable] = buildSetupCatalog(config,setups)
%BUILDSETUPCATALOG Build an N-by-1 car catalog from setup specifications.

if nargin < 2 || isempty(setups)
    setups = config.defaultSetup;
end
validateattributes(config,{'struct'},{'scalar'},mfilename,'config');
validateattributes(setups,{'struct'},{'vector'},mfilename,'setups');
count = numel(setups);
if count == 0
    cars = cell(0,1);
    cases = repmat(emptyCase(),0,1);
    designTable = emptyDesignTable();
    return
end

cars = cell(count,1);
cases = repmat(emptyCase(),count,1);
ids = strings(count,1);
labels = strings(count,1);
sources = strings(count,1);
baselineVersions = strings(count,1);
rearArb = zeros(count,1);
frontSpring = zeros(count,1);
rearSpring = zeros(count,1);
driverWeight = zeros(count,1);
rearWeightDistribution = zeros(count,1);
frontRide = zeros(count,1);
rearRide = zeros(count,1);
aeroMapIds = strings(count,1);
rollSplits = zeros(count,1);
isBaseline = false(count,1);

for i = 1:count
    [cars{i},normalized,derived] = rampSpeed.buildCarFromSetup( ...
        config,setups(i));
    ids(i) = normalized.id;
    labels(i) = normalized.label;
    sources(i) = normalized.source;
    baselineVersions(i) = normalized.baselineVersion;
    rearArb(i) = normalized.rearArbStiffness_NmPerRad;
    frontSpring(i) = normalized.frontSpringRate_lb_in;
    rearSpring(i) = normalized.rearSpringRate_lb_in;
    frontRide(i) = normalized.frontRideHeight_in;
    rearRide(i) = normalized.rearRideHeight_in;
    driverWeight(i) = normalized.driverWeight_kg;
    rearWeightDistribution(i) = normalized.rearWeightDistribution_percent;
    aeroMapIds(i) = normalized.aeroMapId;
    rollSplits(i) = derived.R_sf;
    isBaseline(i) = normalized.isBaseline;

    cases(i).id = normalized.id;
    cases(i).label = normalized.label;
    cases(i).source = normalized.source;
    cases(i).designRow = i;
    cases(i).sourceIndex = i;
    cases(i).carRole = "auto";
    cases(i).carColumn = 1;
    cases(i).setupSpec = normalized;
    cases(i).derived = derived;
    cases(i).isBaseline = isBaseline(i);
end

if numel(unique(ids)) ~= count
    error('rampSpeed:duplicateSetupId', ...
        'Setup IDs must be unique within the ramp-speed catalog.');
end
designTable = table(ids,labels,sources,baselineVersions,rearArb,frontSpring, ...
    rearSpring,frontRide,rearRide,driverWeight,rearWeightDistribution, ...
    aeroMapIds,rollSplits,isBaseline, ...
    'VariableNames',{'id','label','source','baselineVersion', ...
    'rearArbStiffness_NmPerRad','frontSpringRate_lb_in', ...
    'rearSpringRate_lb_in','frontRideHeight_in','rearRideHeight_in', ...
    'driverWeight_kg','rearWeightDistribution_percent', ...
    'aeroMapId','R_sf','isBaseline'});
end

function value = emptyCase()
value = struct('id',"",'label',"",'source',"",'designRow',NaN, ...
    'sourceIndex',NaN,'carRole',"auto",'carColumn',1,'setupSpec',struct(), ...
    'derived',struct(),'isBaseline',false);
end

function value = emptyDesignTable()
value = table(strings(0,1),strings(0,1),strings(0,1),strings(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1), ...
    strings(0,1),zeros(0,1),false(0,1), ...
    'VariableNames',{'id','label','source','baselineVersion', ...
    'rearArbStiffness_NmPerRad','frontSpringRate_lb_in', ...
    'rearSpringRate_lb_in','frontRideHeight_in','rearRideHeight_in', ...
    'driverWeight_kg','rearWeightDistribution_percent', ...
    'aeroMapId','R_sf','isBaseline'});
end
