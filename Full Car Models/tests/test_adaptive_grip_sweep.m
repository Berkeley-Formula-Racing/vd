function tests = test_adaptive_grip_sweep
%TEST_ADAPTIVE_GRIP_SWEEP Contract tests for the adaptive grip entrypoint.
tests = functiontests(localfunctions);
end

function testRejectsMissingTargetEvent(testCase)
targets = struct('autocross',48,'skidpad',4.9);
options = struct('bootstrap',false);

verifyError(testCase,@() run_adaptive_grip_sweep(targets,options), ...
    'run_adaptive_grip_sweep:missingTarget');
end

function testRefinesAroundCandidateSelectsScalarBaselineAndConfirms(testCase)
targets = struct('autocross',48,'skidpad',4.9,'accel',4.8);
calls = cell(0,1);
seenBase = cell(0,1);

options = struct( ...
    'bootstrap',false, ...
    'carConfigFcn',@fixtureCarConfig, ...
    'sweepFcn',@fixtureSweep, ...
    'initialFront',[0.57 0.60 0.63 0.66], ...
    'initialRear',[0.58 0.61 0.64 0.67], ...
    'maxIterations',3, ...
    'minSpan',0.005, ...
    'shrinkFactor',0.5, ...
    'plot',false, ...
    'verbose',false, ...
    'saveResults',false);

result = run_adaptive_grip_sweep(targets,options);

verifyEqual(testCase,seenBase{1},'baseline');
verifyEqual(testCase,numel(result.history),3);
verifyGreaterThan(testCase,result.history(1).frontSpan, ...
    result.history(end).frontSpan);
verifyGreaterThan(testCase,result.history(1).rearSpan, ...
    result.history(end).rearSpan);
verifyEqual(testCase,result.candidate.gripF,0.61,'AbsTol',1e-12);
verifyEqual(testCase,result.candidate.gripR,0.62,'AbsTol',1e-12);
verifyEqual(testCase,result.confirmed.gripF,0.61,'AbsTol',1e-12);
verifyEqual(testCase,result.confirmed.gripR,0.62,'AbsTol',1e-12);
verifyTrue(testCase,numel(calls) > numel(result.history));
verifyEqual(testCase,numel(calls{end}.front),1);
verifyEqual(testCase,numel(calls{end}.rear),1);
    function G = fixtureSweep(baseCell,sweepOptions)
        seenBase{end+1} = baseCell{1,1};
        calls{end+1} = struct('front',sweepOptions.front, ...
            'rear',sweepOptions.rear);

        [F,R] = ndgrid(sweepOptions.front,sweepOptions.rear);
        target = sweepOptions.target;
        autocross = target.autocross + 4*(F-0.61).^2 + 4*(R-0.62).^2;
        skidpad = target.skidpad + 2*(F-0.61).^2 + 5*(R-0.62).^2;
        accel = target.accel + 3*(F-0.61).^2 + 2*(R-0.62).^2;
        score = sqrt(((autocross-target.autocross)/target.autocross/0.01).^2 + ...
                     ((skidpad-target.skidpad)/target.skidpad/0.01).^2 + ...
                     ((accel-target.accel)/target.accel/0.01).^2) / sqrt(3);

        T = table((1:numel(F))',F(:),R(:),skidpad(:),accel(:),autocross(:),score(:), ...
            'VariableNames',{'case','gripF','gripR','skidpad','accel','autocross','score'});
        [~,idx] = min(T.score);

        G = struct();
        G.times = T;
        G.target = target;
        G.settings = struct('front',sweepOptions.front,'rear',sweepOptions.rear, ...
            'events',{{'skidpad','accel','autocross'}});
        G.grid = struct('F',F,'R',R);
        G.best = T(idx,:);
        if numel(sweepOptions.front) > 1 && numel(sweepOptions.rear) > 1
            G.bestInterp = struct('gripF',0.61,'gripR',0.62,'score',0, ...
                'skidpad',target.skidpad,'accel',target.accel, ...
                'autocross',target.autocross,'method','fixture');
        else
            G.bestInterp = struct();
        end
    end
end

function [cars,eventParams,designTable] = fixtureCarConfig()
cars = {'low','lowAccel';'baseline','baselineAccel';'high','highAccel'};
eventParams = struct('fixture',true);
designTable = table([-0.25;0;0.25],[-0.25;0;0.25], ...
    'VariableNames',{'static_front_ride_height_in', ...
                     'static_rear_ride_height_in'});
end
