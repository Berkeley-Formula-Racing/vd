function result = run_lc0_nd_coupling(cfg,makePlots)
%RUN_LC0_ND_COUPLING Fit provisional donor combined-slip coupling behavior.
%   This saves donor-only normalized reference curves and the exponent score;
%   it does not alter target LC0 force levels or vehicle-model code.

if nargin < 1 || isempty(cfg), cfg = lc0NDConfig(); end
if nargin < 2 || isempty(makePlots), makePlots = true; end
required = {'donor','donorLateral','couplingLongitudinalFit', ...
    'couplingLateralFit','couplingFit','outputDirectory'};
if ~all(isfield(cfg,required))
    error('run_lc0_nd_coupling:badConfig', ...
        'CFG requires donor, donorLateral, coupling fit options, and outputDirectory.');
end
if ~isfolder(cfg.outputDirectory), mkdir(cfg.outputDirectory); end

combinedDonor = lc0NDLoadDonor(cfg.donor);
lateralDonor = lc0NDLoadDonor(cfg.donorLateral);
longFit = lc0NDFitLongitudinal(combinedDonor,cfg.couplingLongitudinalFit);
latFit = lc0NDFitLateral(lateralDonor,cfg.couplingLateralFit);
coupling = lc0NDFitCouplingExponent(combinedDonor,longFit,latFit,cfg.couplingFit);
result = struct('config',cfg,'longitudinal',longFit,'lateral',latFit, ...
    'coupling',coupling);
save(fullfile(cfg.outputDirectory,'lc0_nd_donor_coupling.mat'), ...
    'cfg','longFit','latFit','coupling','-v7.3');
if makePlots
    f = figure('Name','Provisional LC0 donor coupling fit','Color','w');
    plot(coupling.score.exponent,coupling.score.rmse_mu,'o-','LineWidth',1.2);
    xlabel('Friction-ellipse exponent'); ylabel('Force coefficient RMSE');
    title(sprintf('Provisional donor coupling fit (%d combined points)', ...
        coupling.n_points));
    grid on;
    exportgraphics(f,fullfile(cfg.outputDirectory,'lc0_nd_coupling_score.png'), ...
        'Resolution',180);
    savefig(f,fullfile(cfg.outputDirectory,'lc0_nd_coupling_score.fig'));
    close(f);
end
end
