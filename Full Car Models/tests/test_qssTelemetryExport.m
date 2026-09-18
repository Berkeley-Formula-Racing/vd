function tests = test_qssTelemetryExport
%TEST_QSSTELEMETRYEXPORT Contract tests for the additive QSS exporter.
tests = functiontests(localfunctions);
end

function testBuildPreservesNativeArraysAndRecordsMapping(testCase)
source = syntheticSource();
result = buildTelemetryResult(source,struct('run_uuid','test-run'));

verifyEqual(testCase,result.manifest.format,'qss-telemetry');
verifyEqual(testCase,result.manifest.schema_version,'1.0');
verifyEqual(testCase,result.cases(1).laps(1).channels(1).values(:), ...
    [10;11;12;13;14],'AbsTol',0);
verifyEqual(testCase,result.cases(1).laps(1).channels(2).values(:), ...
    [0.1;0.2;0.3],'AbsTol',0);
verifyTrue(testCase,isfield(result.cases(1).laps(1).metadata, ...
    'node_segment_mapping'));
verifyEqual(testCase,result.cases(1).laps(1).metadata.time_semantics, ...
    'raw_native_time');
verifyTrue(testCase,validateQSSResults(result));
end

function testExportWritesExpectedPathsAndUtf8Metadata(testCase)
source = syntheticSource();
result = buildTelemetryResult(source,struct('run_uuid','utf8-run'));
folder = tempname;
mkdir(folder);
cleanup = onCleanup(@() rmdir(folder,'s')); %#ok<NASGU>
path = fullfile(folder,'telemetry.h5');
actual = exportQSSResults(result,path);

verifyEqual(testCase,string(actual),string(path));
verifyEqual(testCase,h5read(actual,'/manifest_json'), ...
    uint8(unicode2native(jsonencode(result.manifest),'UTF-8'))');
verifyEqual(testCase,h5read(actual, ...
    '/cases/baseline/laps/flying/channels/speed_mps/values'), ...
    [10;11;12;13;14],'AbsTol',0);
verifyEqual(testCase,h5read(actual, ...
    '/cases/baseline/laps/flying/channels/steering_angle_rad/values'), ...
    [0.1;0.2;0.3],'AbsTol',0);
metadata = native2unicode(h5read(actual, ...
    '/cases/baseline/metadata_json'),'UTF-8');
beta = native2unicode(uint8([206 178]),'UTF-8');
verifyTrue(testCase,any(contains(string(metadata),beta)));
end

function testReconstructionFlagsUnavailablePointsWithoutInterpolation(testCase)
source = syntheticSource();
result = buildTelemetryResult(source,struct('reconstruct_detail',true));
lap = result.cases(1).laps(1);
detail = reconstructLapTelemetry([],lap,struct('spacing_m',2));

verifyTrue(testCase,all(strcmp(detail.diagnostics.status,'unavailable')));
verifyTrue(testCase,all(~detail.diagnostics.valid));
verifyTrue(testCase,all(isnan(detail.channels(1).values)));
verifyEqual(testCase,detail.diagnostics.failure_reason{1},'car_model_unavailable');
end

function testValidatorRejectsFiniteInvalidSamples(testCase)
source = syntheticSource();
result = buildTelemetryResult(source);
result.cases(1).laps(1).channels(1).valid(2) = false;
result.cases(1).laps(1).channels(1).values(2) = 42;

[valid,report] = validateQSSResults(result);
verifyFalse(testCase,valid);
verifyTrue(testCase,contains(string(report.errors{1}),'invalid'));
end

function source = syntheticSource()
source = struct();
source.source_type = 'plain_study';
beta = native2unicode(uint8([206 178]),'UTF-8');
source.track = struct('id','unit_track','metadata',struct( ...
    'name',[beta 'eta track'],'geometry_source','curvature'), ...
    'distance_m',[0;2;4;6;8], ...
    'curvature_per_m',[0;0.02;0.02;0;-0.01]);
source.case_id = 'baseline';
source.case_metadata = struct('label',[beta 'eta case']);
source.setup = struct('mass_kg',250);
source.laps = struct();
source.laps(1).id = 'flying';
source.laps(1).event = 'autocross';
source.laps(1).role = 'flying';
source.laps(1).track_id = 'unit_track';
source.laps(1).time_s = [0;1;2;3;4];
source.laps(1).distance_m = [0;2;4;6;8];
source.laps(1).channels = struct(...
    'id',{'speed_mps','steering_angle_rad'}, ...
    'label',{'Vehicle speed','Steering'}, ...
    'unit',{'m/s','rad'}, ...
    'values',{[10;11;12;13;14],[0.1;0.2;0.3]}, ...
    'origin',{'simulation_output','simulation_output'}, ...
    'interpolation',{'linear','linear'}, ...
    'description',{'',''}, ...
    'coordinate_frame',{'vehicle','vehicle'}, ...
    'sign_convention',{'forward positive','left positive'});
end
