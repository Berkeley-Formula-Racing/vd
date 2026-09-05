function result = run_lc0_nd_prototype(cfg,makePlots)
%RUN_LC0_ND_PROTOTYPE Summarise measured LC0 lateral data for calibration.
%   RESULT = RUN_LC0_ND_PROTOTYPE() loads the configured free-rolling LC0
%   target runs, writes their SI-unit summary to results, and saves two
%   diagnostic figures.  It does not alter Tire2, carConfig, or any fitted
%   vehicle-model parameter.

if nargin < 1 || isempty(cfg), cfg = lc0NDConfig(); end
if nargin < 2 || isempty(makePlots), makePlots = true; end
if ~isfield(cfg,'summary')
    cfg.summary = defaultSummaryOptions();
end
if ~isfolder(cfg.outputDirectory)
    mkdir(cfg.outputDirectory);
end

[data,manifest] = lc0NDLoadFreeRolling(cfg);
summary = lc0NDLateralSummary(data,cfg.summary);
result = struct('config',cfg,'manifest',manifest,'data',data,'summary',summary);
save(fullfile(cfg.outputDirectory,'lc0_nd_target_summary.mat'), ...
    'cfg','manifest','data','summary','-v7.3');

if makePlots
    saveDiagnosticPlots(result,cfg.outputDirectory);
end
end

function opts = defaultSummaryOptions()
lbfToN = 4.4482216152605;
opts = struct('loadCenters_N',(75:50:275)'*lbfToN, ...
    'loadHalfWidth_N',15*lbfToN, ...
    'pressureCenters_psi',[8;10;12;14], ...
    'pressureHalfWidth_psi',1.0, ...
    'camberCenters_deg',[0;2;4], ...
    'camberHalfWidth_deg',0.6, ...
    'smallSlipWindow_deg',1, ...
    'minSamples',100, ...
    'minSlipAngleSpan_deg',8);
end

function saveDiagnosticPlots(result,outputDirectory)
summary = result.summary;
valid = summary.is_complete_sweep & isfinite(summary.peak_mu_y);
if any(valid)
    f = figure('Name','LC0 measured lateral summary','Color','w');
    scatter(summary.mean_Fz_N(valid),summary.peak_mu_y(valid),55, ...
        summary.mean_pressure_psi(valid),'filled');
    xlabel('Normal load [N]'); ylabel('Measured peak lateral \mu');
    title('16x7.5-10 LC0 free-rolling lateral measurements');
    grid on; cb = colorbar; cb.Label.String = 'Pressure [psi]';
    exportgraphics(f,fullfile(outputDirectory,'lc0_nd_peak_mu_vs_load.png'), ...
        'Resolution',180);
    savefig(f,fullfile(outputDirectory,'lc0_nd_peak_mu_vs_load.fig'));
    close(f);

    f = figure('Name','LC0 measured lateral stiffness','Color','w');
    scatter(summary.mean_Fz_N(valid),summary.stiffness_N_per_deg(valid),55, ...
        summary.mean_camber_deg(valid),'filled');
    xlabel('Normal load [N]'); ylabel('Small-angle stiffness [N/deg]');
    title('16x7.5-10 LC0 free-rolling lateral stiffness');
    grid on; cb = colorbar; cb.Label.String = 'Camber [deg]';
    exportgraphics(f,fullfile(outputDirectory,'lc0_nd_stiffness_vs_load.png'), ...
        'Resolution',180);
    savefig(f,fullfile(outputDirectory,'lc0_nd_stiffness_vs_load.fig'));
    close(f);
end

f = figure('Name','LC0 data coverage','Color','w');
scatter(result.data.Fz_N,result.data.pressure_psi,8,result.data.camber_deg,'filled');
xlabel('Normal load [N]'); ylabel('Pressure [psi]');
title('LC0 free-rolling data coverage'); grid on;
cb = colorbar; cb.Label.String = 'Camber [deg]';
exportgraphics(f,fullfile(outputDirectory,'lc0_nd_data_coverage.png'), ...
    'Resolution',180);
savefig(f,fullfile(outputDirectory,'lc0_nd_data_coverage.fig'));
close(f);
end
