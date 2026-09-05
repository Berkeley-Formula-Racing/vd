function result = run_lc0_nd_donor(cfg,makePlots)
%RUN_LC0_ND_DONOR Characterise the provisional LC0 donor longitudinal data.
%   The saved curve is normalized by donor normal load and remains separate
%   from the target-car tyre model until vehicle-data validation is done.

if nargin < 1 || isempty(cfg), cfg = lc0NDConfig(); end
if nargin < 2 || isempty(makePlots), makePlots = true; end
if ~isfield(cfg,'donor') || ~isfield(cfg,'donorFit') || ...
        ~isfield(cfg,'outputDirectory')
    error('run_lc0_nd_donor:badConfig', ...
        'CFG requires donor, donorFit, and outputDirectory fields.');
end
if ~isfolder(cfg.outputDirectory)
    mkdir(cfg.outputDirectory);
end

donor = lc0NDLoadDonor(cfg.donor);
fit = lc0NDFitLongitudinal(donor,cfg.donorFit);
result = struct('config',cfg,'donor',donor,'fit',fit);
save(fullfile(cfg.outputDirectory,'lc0_nd_donor_longitudinal.mat'), ...
    'cfg','donor','fit','-v7.3');
if makePlots
    plotDonorLongitudinal(fit,cfg.outputDirectory);
end
end

function plotDonorLongitudinal(fit,outputDirectory)
curve = fit.curve;
good = curve.is_qualified;
f = figure('Name','Provisional LC0 donor longitudinal reference','Color','w');
errorbar(curve.slip_ratio(good),curve.mu_x(good),curve.mu_x_std(good), ...
    'o-','LineWidth',1.2,'MarkerFaceColor',[0.1 0.45 0.8]);
xlabel('Slip ratio'); ylabel('Donor longitudinal coefficient, \mu_x');
title(sprintf('Provisional donor longitudinal reference (%d selected samples)', ...
    fit.n_selected));
grid on;
exportgraphics(f,fullfile(outputDirectory,'lc0_nd_donor_longitudinal.png'), ...
    'Resolution',180);
savefig(f,fullfile(outputDirectory,'lc0_nd_donor_longitudinal.fig'));
close(f);
end
