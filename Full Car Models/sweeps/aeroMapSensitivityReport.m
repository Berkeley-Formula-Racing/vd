function aeroMapSensitivityReport(S)
%AEROMAPSENSITIVITYREPORT Print sensitivity and aeromap operating diagnostics.

fprintf('\n=====================================================================\n');
fprintf(' AEROMAP SENSITIVITY SUMMARY\n');
fprintf('=====================================================================\n');
fprintf('\nCases (baseline marked *):\n');
fprintf('  %-5s %-20s %-10s %9s %9s\n','case','input','value','F RH (in)','R RH (in)');
for i = 1:height(S.cases)
    mark = ' '; if S.cases.isBaseline(i), mark = '*'; end
    fprintf('%s %-5d %-20s %+10.4f %9.4f %9.4f\n',mark,S.cases.case(i), ...
        S.cases.swept_parameter(i),S.cases.swept_value(i), ...
        S.cases.front_static_ride_height_in(i),S.cases.rear_static_ride_height_in(i));
end

fprintf('\nSolve and aeromap coverage:\n');
for i = 1:height(S.cases)
    rows = S.metrics(S.metrics.case == S.cases.case(i),:);
    converged = mean(rows.exitflag == 1 | rows.exitflag == 2)*100;
    outside = mean(rows.aero_outside_map)*100;
    residual = max(rows.aero_residual_in,[],'omitnan');
    fprintf('  case %-3d %6.1f%% converged, %5.1f%% outside map, max residual %.2g in\n', ...
        S.cases.case(i),converged,outside,residual);
end

for i = 1:numel(S.paramInfo)
    spec = S.paramInfo(i);
    P = S.sens.(char(spec.name));
    fprintf('\n---------------------------------------------------------------------\n');
    fprintf(' d(metric)/d%s   baseline = %.4f %s\n',spec.label,P.pBase,spec.unit);
    fprintf(' perturbations: %s\n',mat2str(round(P.pValues(:)',4)));
    fprintf('  peak capability, averaged over velocity:\n');
    reportMean(P,'maxGLat_dp','max lateral g');
    reportMean(P,'maxGLong_dp','max acceleration g');
    reportMean(P,'maxGBrake_dp','max braking g');
    if isfield(P,'d_autocross_dp')
        fprintf('  autocross: d(time)/d%s = %+9.4f s per %s\n', ...
            spec.label,P.d_autocross_dp,spec.unit);
    end
end
fprintf('\n=====================================================================\n');
end

function reportMean(P,field,label)
values = P.(field);
fprintf('    %-22s %+8.4f  (range %+.4f to %+.4f)\n',label, ...
    mean(values,'omitnan'),min(values,[],'omitnan'),max(values,[],'omitnan'));
end
