load("C:\Users\johny\Documents\VD\600DOE.mat");

%% --- 1. DEFINE THE DYNAMIC VARIABLE MAP ---
% Format: {'RSM_Name', 'Car_Object_Property_Path'}
% This map includes everything from your carConfig except engine/gears.
varMap = {
    % Chassis
    'mass_total',   'M';
    'wheelbase',    'W_b';
    'cg_height',    'h_g'; 
    'track_width',  't_f';
    'roll_stiff_f', 'R_sf';
    'izz',          'I_zz';
    'ackermann',    'ackermann';
    'compliance',   'camber_compliance';
    
    % Aero (Main Car) - UPDATED TO MATCH AERO CLASS PROPERTIES
    'cla',          'aero.cla';  % Double-check if this is ClA or ClA_tot
    'cda',          'aero.cda';  % Double-check if this is CdA or CdA_tot
    'aero_dist',    'aero.D_f';  % Changed from .distribution to .D_f
    
    % Drivetrain & Brakes
    'final_drive',  'powertrain.final_drive';
    'drive_eff',    'powertrain.drivetrain_efficiency';
    'tbr_gain',     'powertrain.G_d2_driving';
    'brake_dist',   'powertrain.brake_distribution';
    'max_brake',    'powertrain.max_braking_torque';
    
    % Tires
    'static_gamma', 'static_gamma';
    'pressure',     'tire.p_i';
    'mu_scale',     'tire.friction_scaling_factor'
};

% Special case: Variables that only exist in the Accel Car (Column 2)
accelVarMap = {
    'accel_cda',    'aero.cda';
    'accel_cla',    'aero.cla'
};

%% --- 2. DYNAMIC EXTRACTION LOOP ---
numRuns = size(carCell, 1);
outputNames = {'t_autox', 't_accel', 't_skid', 'total_work'};

% Initialize data storage
inputData = zeros(numRuns, size(varMap, 1) + size(accelVarMap, 1));
outputData = NaN(numRuns, numel(outputNames));

for i = 1:numRuns
    carMain = carCell{i, 1};
    carAccel = carCell{i, 2};
    
    % --- Extract Main Car Variables ---
    for v = 1:size(varMap, 1)
        inputData(i, v) = getProp(carMain, varMap{v, 2});
    end
    
    % --- Extract Accel-Specific Variables ---
    offset = size(varMap, 1);
    for v = 1:size(accelVarMap, 1)
        inputData(i, offset + v) = getProp(carAccel, accelVarMap{v, 2});
    end
    
    % --- Extract Performance Results ---
    if ~isempty(carMain.comp)
        c = carMain.comp;
        outputData(i, :) = [c.times.autocross, c.times.accel, c.times.skidpad, ...
                            c.autocross.metrics.total_work_J];
    end
end

% Combine into a Table
allNames = [varMap(:, 1)', accelVarMap(:, 1)', outputNames];
doeTable = array2table([inputData, outputData], 'VariableNames', allNames);

% Clean: Remove failed runs and variables with zero variance (constants)
doeTable = doeTable(~isnan(doeTable.t_autox), :);
doeTable = doeTable(:, var(table2array(doeTable)) > 1e-10 | ismember(doeTable.Properties.VariableNames, outputNames));


responseVars = {'t_autox', 't_accel', 't_skid', 'total_work'};

for i = 1:numel(responseVars)
    varName = responseVars{i};
    % Force the column to be real numbers only
    doeTable.(varName) = real(doeTable.(varName));
end
% --- DATA SANITIZATION & GHOST FILTERING ---
% 1. Strip Numerical Noise (Complex Numbers)
% StepwiseLM and bar charts will crash if even 1e-18 of an imaginary component exists.
responseVars = {'t_autox', 't_accel', 't_skid', 'total_work'};
for i = 1:numel(responseVars)
    doeTable.(responseVars{i}) = real(doeTable.(responseVars{i}));
end

% 2. Define "Physicality Window"
% Based on Michigan-style track (1km) and standard 8.5m Skidpad
limits.skid_max  = 8.0;   % Anything over 7s is a solver stall/crawl
limits.autox_max = 75.0;  % Upper bound for a "slow" but valid lap
limits.autox_min = 35.0;  % Catch "teleporting" cars (10s laps are impossible)
limits.accel_min = 3.0;   % Catch failed launches

% 3. Identify Logic Failures
% We create separate masks to see exactly WHY cars are failing
is_too_slow_skid  = doeTable.t_skid  > limits.skid_max;
is_too_fast_autox = doeTable.t_autox < limits.autox_min;
is_too_slow_autox = doeTable.t_autox > limits.autox_max;
is_broken_accel   = doeTable.t_accel < limits.accel_min;

% Combine into a master "Ghost" mask
is_ghost = is_too_slow_skid | is_too_fast_autox | is_too_slow_autox | is_broken_accel;

% 4. Logging & Diagnostics
fprintf('\n--- DOE Data Quality Report ---\n');
if any(is_ghost)
    fprintf('  > Detected %d "Ghost" setups:\n', sum(is_ghost));
    if sum(is_too_fast_autox) > 0, fprintf('    - Autocross Teleports (<35s): %d\n', sum(is_too_fast_autox)); end
    if sum(is_too_slow_skid)  > 0, fprintf('    - Skidpad Stalls (>8s):      %d\n', sum(is_too_slow_skid)); end
    if sum(is_too_slow_autox) > 0, fprintf('    - Autocross DNFs (>75s):     %d\n', sum(is_too_slow_autox)); end
else
    fprintf('  > All runs passed physicality checks.\n');
end

% 5. Apply Filter
doeTable = doeTable(~is_ghost, :);

fprintf('  > Final Result: %d valid runs remaining for RSM fitting.\n', height(doeTable));
fprintf('  > Variables Tracked: %d\n', width(doeTable) - numel(responseVars));
fprintf('-------------------------------\n');

%% --- 3. AUTOMATED STEPWISE RSM FITTING ---
targets = {'t_autox', 't_accel', 't_skid', 'total_work'};
% Base inputs (everything except the target lap times)
inputVars = setdiff(doeTable.Properties.VariableNames, outputNames, 'stable');

% --- MANUAL OVERRIDE: EXCLUDE VARIABLES FROM RSM ---
% Add any variable names here that you want the model to ignore.
varsToExclude = {}; 
inputVars = setdiff(inputVars, varsToExclude, 'stable');
models = struct();

for t = 1:numel(targets)
    resp = targets{t};
    fprintf('\n--- Analyzing %s ---\n', resp);
    
    % Filter table for target + relevant inputs
    currTable = doeTable(:, [inputVars, {resp}]);
    
    % STEPWISE: Automatically picks only the parameters that matter
    models.(resp) = stepwiselm(currTable, 'ResponseVar', resp, ...
        'Lower', 'constant', 'Upper', 'quadratic', 'Criterion', 'bic');
    
    fprintf('Significant Predictors found: %d\n', numel(models.(resp).PredictorNames));
    fprintf('Adjusted R-squared: %.4f\n', models.(resp).Rsquared.Adjusted);
end

%% --- 4. VISUALIZATION ---
figure('Name', 'Design Sensitivity', 'Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.7]);
subplot(1,4,1); plotSlice(models.t_autox); title('Autocross Sensitivity');
subplot(1,4,2); plotSlice(models.t_accel); title('Accel Sensitivity');
subplot(1,4,3); plotSlice(models.t_skid);  title('Skidpad Sensitivity');
subplot(1,4,4); plotSlice(models.total_work);  title('Total Work Sensitivity');
%% --- HELPER FUNCTION: Recursive Property Search ---
function val = getProp(obj, pathStr)
    path = strsplit(pathStr, '.');
    val = obj;
    for p = 1:numel(path)
        val = val.(path{p});
    end
end

%% --- EVENT SENSITIVITY VISUALIZATION ---
events = {'t_autox', 't_accel', 't_skid', 'total_work'};
colors = {[0.2 0.4 0.8], [0.8 0.2 0.2], [0.2 0.4 0.8], [0.2 0.7 0.2]}; % Blue, Red, Green

for e = 1:numel(events)
    target = events{e};
    mdl = models.(target);
    
    % 1. CHECK FOR EMPTY MODEL
    % If BIC kicked everything out, names will be empty.
    if numel(mdl.CoefficientNames) <= 1
        fprintf('Warning: No significant factors found for %s. Skipping plot.\n', target);
        continue; 
    end
    
    meanVal = mean(doeTable.(target));
    coeffs = mdl.Coefficients.Estimate(2:end);
    names = mdl.CoefficientNames(2:end);
    
    % 2. Calculate the "Swing" (Impact) of each term
    % We want to know: (Beta * Range) / MeanTime * 100
    percentImpact = zeros(numel(coeffs), 1);
    
    for c = 1:numel(coeffs)
        termName = names{c};
        
        % Check if it's an interaction (e.g., 'mass:cla')
        if contains(termName, ':')
            parts = strsplit(termName, ':');
            % Range impact of interaction: Beta * (maxA*maxB - minA*minB)
            rangeA = [min(doeTable.(parts{1})), max(doeTable.(parts{1}))];
            rangeB = [min(doeTable.(parts{2})), max(doeTable.(parts{2}))];
            valDelta = (rangeA(2)*rangeB(2)) - (rangeA(1)*rangeB(1));
        else
            % Linear or Quadratic term (e.g., 'mass' or 'mass^2')
            % Stripping the '^2' for range calculation
            cleanName = regexprep(termName, '\^2', '');
            valDelta = max(doeTable.(cleanName)) - min(doeTable.(cleanName));
            if contains(termName, '^2')
                valDelta = max(doeTable.(cleanName))^2 - min(doeTable.(cleanName))^2;
            end
        end
        
        % Calculate Percent Change
        percentImpact(c) = (coeffs(c) * valDelta / meanVal) * 100;
    end
    
    % 2. SORT AND PLOT
    [~, sortIdx] = sort(abs(percentImpact), 'descend');
    sortedImpact = percentImpact(sortIdx);
    sortedNames = names(sortIdx);
    
    figure('Name', ['Sensitivity: ' target], 'Color', 'w');
    h = barh(sortedImpact);
    h.FaceColor = colors{e};
    
    % 3. THE FIX: Explicitly set ticks and turn off the TeX interpreter
    % This prevents "final_drive" from trying to render as "final_{drive}"
    ax = gca;
    set(ax, 'YTick', 1:numel(sortedNames));
    set(ax, 'YTickLabel', sortedNames);
    set(ax, 'TickLabelInterpreter', 'none'); % CRITICAL: Disables auto-subscripting
    set(ax, 'YDir', 'reverse');
    
    xlabel('Total Effect on Metric (%)');
    title(['Primary Drivers: ', strrep(target, '_', ' ')]);
    grid on;
    
    % 4. BAR LABELS
    for i = 1:numel(sortedImpact)
        text(sortedImpact(i), i, [' ', num2str(sortedImpact(i), '%+0.2f'), '%'], ...
            'VerticalAlignment', 'middle', 'Interpreter', 'none');
    end
end

%% --- 5. 3D COMPOUND EFFECTS (INTERACTION) VISUALIZATION ---
% Choose your target event and the two factors you want to investigate
targetEvent = 't_autox';      % Options: 't_autox', 't_accel', 't_skid', 'total_work'
factor1     = 'mass_total';   % Must match a variable name in your doeTable
factor2     = 'cla';          % Must match a variable name in your doeTable

% Grab the specific Stepwise Model and original table
mdl = models.(targetEvent);
predictors = mdl.PredictorNames;

% Ensure the factors were actually kept by the stepwise model
if ~ismember(factor1, predictors) || ~ismember(factor2, predictors)
    warning('One or both factors were excluded by the Stepwise model because they lacked significance for %s.', targetEvent);
end

% 1. Setup the Evaluation Grid (50x50 resolution)
% Find the min and max limits from your DOE table
f1_min = min(doeTable.(factor1)); f1_max = max(doeTable.(factor1));
f2_min = min(doeTable.(factor2)); f2_max = max(doeTable.(factor2));

[X1_grid, X2_grid] = meshgrid(linspace(f1_min, f1_max, 50), ...
                              linspace(f2_min, f2_max, 50));

% 2. Create "Dummy Data" holding all other variables at their mean
% This isolates the effects to JUST the two variables we care about
numPoints = numel(X1_grid);
dummyData = table();
for i = 1:numel(predictors)
    varName = predictors{i};
    % Fill with the mean value across all points
    meanVal = mean(doeTable.(varName));
    dummyData.(varName) = repmat(meanVal, numPoints, 1);
end

% 3. Inject our sweeping 3D grid variables into the table
dummyData.(factor1) = X1_grid(:);
dummyData.(factor2) = X2_grid(:);

% 4. Predict the response using the fitted RSM
Z_pred = predict(mdl, dummyData);
Z_grid = reshape(Z_pred, size(X1_grid));

% 5. Render the 3D Surface Plot
figure('Name', sprintf('Interaction: %s vs %s', factor1, factor2), 'Color', 'w');
surf(X1_grid, X2_grid, Z_grid, 'EdgeColor', 'none', 'FaceAlpha', 0.85);
colormap parula;
colorbar;

% Plot Formatting
ax = gca;
set(ax, 'TickLabelInterpreter', 'none'); 
xlabel(factor1, 'Interpreter', 'none', 'FontWeight', 'bold');
ylabel(factor2, 'Interpreter', 'none', 'FontWeight', 'bold');
zlabel(strrep(targetEvent, '_', ' '), 'Interpreter', 'none', 'FontWeight', 'bold');
title(sprintf('Compound Effect on %s\n(Other variables held at mean)', strrep(targetEvent, '_', ' ')));

% Adjust view angle for standard 3D perspective
view(-45, 30);
grid on;