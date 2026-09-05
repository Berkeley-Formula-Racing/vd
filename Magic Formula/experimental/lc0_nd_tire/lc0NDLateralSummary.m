function summary = lc0NDLateralSummary(data,opts)
%LC0NDLATERALSUMMARY Summarise free-rolling LC0 lateral measurements by bin.
%   SUMMARY = LC0NDLATERALSUMMARY(DATA,OPTS) reports measured peak lateral
%   force, peak friction and small-angle stiffness for every requested
%   load/pressure/camber bin.  DATA uses the positive-normal-load SI sign
%   convention produced by lc0NDLoadFreeRolling.

requiredData = {'slipAngle_deg','Fy_N','Fz_N','pressure_psi','camber_deg'};
for i = 1:numel(requiredData)
    if ~isfield(data,requiredData{i})
        error('lc0NDLateralSummary:missingField', ...
            'DATA must contain %s.',requiredData{i});
    end
end
requiredOpts = {'loadCenters_N','loadHalfWidth_N','pressureCenters_psi', ...
    'pressureHalfWidth_psi','camberCenters_deg','camberHalfWidth_deg', ...
    'smallSlipWindow_deg'};
for i = 1:numel(requiredOpts)
    if ~isfield(opts,requiredOpts{i})
        error('lc0NDLateralSummary:missingOption', ...
            'OPTS must contain %s.',requiredOpts{i});
    end
end

sa = data.slipAngle_deg(:);
fy = data.Fy_N(:);
fz = data.Fz_N(:);
pressure = data.pressure_psi(:);
camber = data.camber_deg(:);
n = numel(sa);
if any([numel(fy),numel(fz),numel(pressure),numel(camber)] ~= n)
    error('lc0NDLateralSummary:sizeMismatch', ...
        'All DATA channels must have equal length.');
end
minSamples = getOption(opts,'minSamples',1);
minSlipAngleSpan = getOption(opts,'minSlipAngleSpan_deg',0);

loadCenters = opts.loadCenters_N(:);
pressureCenters = opts.pressureCenters_psi(:);
camberCenters = opts.camberCenters_deg(:);
nRows = numel(loadCenters)*numel(pressureCenters)*numel(camberCenters);
rows = nan(nRows,13);
row = 0;
for iLoad = 1:numel(loadCenters)
    for iPressure = 1:numel(pressureCenters)
        for iCamber = 1:numel(camberCenters)
            row = row + 1;
            use = isfinite(sa) & isfinite(fy) & isfinite(fz) & ...
                isfinite(pressure) & isfinite(camber) & ...
                abs(fz-loadCenters(iLoad)) <= opts.loadHalfWidth_N & ...
                abs(pressure-pressureCenters(iPressure)) <= opts.pressureHalfWidth_psi & ...
                abs(camber-camberCenters(iCamber)) <= opts.camberHalfWidth_deg;
            rows(row,1:3) = [loadCenters(iLoad),pressureCenters(iPressure), ...
                camberCenters(iCamber)];
            rows(row,4) = sum(use);
            if ~any(use)
                continue
            end

            saBin = sa(use);
            fyBin = fy(use);
            fzBin = fz(use);
            pressureBin = pressure(use);
            camberBin = camber(use);
            [peakAbs,peakIndex] = max(abs(fyBin));
            rows(row,5:10) = [mean(fzBin),mean(pressureBin),mean(camberBin), ...
                peakAbs,peakAbs/mean(fzBin),saBin(peakIndex)];
            smallSlip = abs(saBin) <= opts.smallSlipWindow_deg;
            denom = sum(saBin(smallSlip).^2);
            if denom > 0
                rows(row,11) = abs(sum(saBin(smallSlip).*fyBin(smallSlip))/denom);
            end
            rows(row,12) = max(saBin)-min(saBin);
            rows(row,13) = rows(row,4) >= minSamples && ...
                rows(row,12) >= minSlipAngleSpan;
        end
    end
end

summary = array2table(rows,'VariableNames', ...
    {'load_center_N','pressure_center_psi','camber_center_deg','n_samples', ...
    'mean_Fz_N','mean_pressure_psi','mean_camber_deg','peak_abs_Fy_N', ...
    'peak_mu_y','alpha_at_peak_deg','stiffness_N_per_deg', ...
    'slipAngle_span_deg','is_complete_sweep'});
complete = summary.is_complete_sweep;
complete(~isfinite(complete)) = 0;
summary.is_complete_sweep = logical(complete);
end

function value = getOption(opts,name,defaultValue)
if isfield(opts,name)
    value = opts.(name);
else
    value = defaultValue;
end
if ~isscalar(value) || ~isfinite(value) || value < 0
    error('lc0NDLateralSummary:badOption', ...
        'OPTS.%s must be a finite, nonnegative scalar.',name);
end
end
