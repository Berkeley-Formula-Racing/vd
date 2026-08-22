function aeroSensitivityReport(S)
% Prints a text summary of an aeroSensitivity result, including the health of
% the underlying solves. Read the convergence line first: sensitivities built
% on cases that did not converge are not worth interpreting.

fprintf('\n=====================================================================\n');
fprintf(' AERO SENSITIVITY SUMMARY\n');
fprintf('=====================================================================\n');

fprintf('\nCases (baseline marked *):\n');
fprintf('  %-5s %-9s %-9s %-9s\n','case','ClA','CdA','CoP');
for i = 1:height(S.cases)
    mark = ' '; if S.cases.isBaseline(i), mark = '*'; end
    fprintf('%s %-5d %-9.4f %-9.4f %-9.4f\n',mark,S.cases.case(i), ...
        S.cases.ClA(i),S.cases.CdA(i),S.cases.CoP(i));
end

% solve health
fprintf('\nSolve convergence (share of g-g points with exitflag 1 or 2):\n');
for i = 1:height(S.cases)
    m = S.metrics(S.metrics.case==i,:);
    if isempty(m), fprintf('  case %-3d no points\n',i); continue, end
    ok = mean(m.exitflag==1 | m.exitflag==2)*100;
    lift = mean(m.min_Fz < 0)*100;
    fprintf('  case %-3d %6.1f%% converged, %5.1f%% of points at wheel lift, %d points\n', ...
        i,ok,lift,height(m));
end

params = fieldnames(S.sens);
for p = 1:numel(params)
    name = params{p};
    P = S.sens.(name);
    fprintf('\n---------------------------------------------------------------------\n');
    fprintf(' d(metric)/d%s     baseline %s = %.4f\n',name,name,P.pBase);
    fprintf(' perturbation levels: %s\n',mat2str(round(P.pValues(:)',4)));
    fprintf('%s\n',P.note);
    fprintf('---------------------------------------------------------------------\n');

    fprintf('  peak capability, averaged over velocity:\n');
    rep = {'maxGLat_dp','max lateral g'; ...
           'maxGLong_dp','max acceleration g'; ...
           'maxGBrake_dp','max braking g'};
    for r = 1:size(rep,1)
        y = P.(rep{r,1});
        fprintf('    %-22s %+8.4f  (range over v: %+.4f to %+.4f)\n', ...
            rep{r,2},mean(y,'omitnan'),min(y),max(y));
    end

    % nonlinearity check: do the perturbation levels agree on the slope?
    sp = P.maxGLat_dp_spread;
    sl = abs(P.maxGLat_dp);
    rel = mean(sp(isfinite(sp)))/max(mean(sl(isfinite(sl))),eps);
    if rel > 0.25
        fprintf('    NOTE: slope varies %.0f%% across perturbation levels --\n',rel*100);
        fprintf('          response is nonlinear, do not read a single gradient.\n');
    end

    % reliability of the theta-resolved envelope sensitivity
    solved = isfinite(P.r_theta_dp);
    if any(solved(:))
        bad = solved & ~P.r_theta_reliable;
        fprintf('  envelope vs theta: %.0f%% of solved cells reliable',...
            100*(1-sum(bad(:))/sum(solved(:))));
        if any(bad(:))
            vIdx = any(bad,2);
            fprintf(' (sign disagreement at v = %s m/s)', ...
                mat2str(round(reshape(S.grids.vCar(vIdx),1,[]),1)));
        end
        fprintf('\n');
        fprintf('    unreliable cells are blanked in the continuous plots\n');
    end

    tf = fieldnames(P);
    tf = tf(startsWith(tf,'d_') & endsWith(tf,'_dp'));
    if ~isempty(tf)
        fprintf('  lap time (negative = faster):\n');
        for t = 1:numel(tf)
            % anchored strip: erase(...,'d_') would also eat the 'd_' inside
            % names like skidpad_dp
            ev = regexprep(tf{t},'^d_(.*)_dp$','$1');
            fprintf('    d(%-10s)/d%-4s = %+9.4f s per unit\n',ev,name,P.(tf{t}));
        end
    end
end
fprintf('\n=====================================================================\n');
end
