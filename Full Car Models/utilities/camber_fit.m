%% Camber response functions using separate wheel steer inputs
clear; clc;

% --- 1) Load & rename -----------------------------------------------------
t = readtable("C:\Users\johny\Documents\VD\vd\winGeo\front_rs.CSV");
t = renamevars(t, ["Var3","Var4","Var5","Var6","Var7","Var8"], ...
                  ["roll","steer","CamberL","CamberR","SteerL","SteerR"]);

t = t(:, {'roll','CamberL','CamberR','SteerL','SteerR'});


mask = all(isfinite(t{:,:}),2);
t = t(mask,:);

% --- 2) Quadratic design terms for each side ------------------------------
% Left
tblL = table(t.roll, t.SteerL, t.roll.^2, t.SteerL.^2, t.roll.*t.SteerL, t.CamberL, ...
             'VariableNames', {'roll','steer','roll2','steer2','cross','Camber'});
% Right
tblR = table(t.roll, t.SteerR, t.roll.^2, t.SteerR.^2, t.roll.*t.SteerR, t.CamberR, ...
             'VariableNames', {'roll','steer','roll2','steer2','cross','Camber'});

% --- 3) Fit quadratic models ----------------------------------------------
mdlL = fitlm(tblL, 'Camber ~ roll + steer + roll2 + steer2 + cross', ...
                    'Intercept', true, 'RobustOpts','on');
mdlR = fitlm(tblR, 'Camber ~ roll + steer + roll2 + steer2 + cross', ...
                    'Intercept', true, 'RobustOpts','on');

disp('=== Left camber model ===');  disp(mdlL);
disp('=== Right camber model ==='); disp(mdlR);

% --- 4) Function handles --------------------------------------------------
camberL_fun = @(roll,steerL) local_predict(mdlL, roll, steerL);
camberR_fun = @(roll,steerR) local_predict(mdlR, roll, steerR);

% --- 5) Quick check at data points ----------------------------------------
yL_hat = camberL_fun(t.roll, t.SteerL);
yR_hat = camberR_fun(t.roll, t.SteerR);

rmseL = sqrt(mean((yL_hat - t.CamberL).^2));
rmseR = sqrt(mean((yR_hat - t.CamberR).^2));

%% -------------------- local function -------------------------------------
function y = local_predict(mdl, roll, steer)
    sz = size(roll + steer);
    rr = roll(:); ss = steer(:);
    TT = table(rr, ss, rr.^2, ss.^2, rr.*ss, ...
               'VariableNames', {'roll','steer','roll2','steer2','cross'});
    y = predict(mdl, TT);
    y = reshape(y, sz);
end
% After you fit mdlL / mdlR:
getb = @(m,term) m.Coefficients.Estimate(strcmp(m.CoefficientNames,term));

betaL = [ ...
    getb(mdlL,'(Intercept)'); ...
    getb(mdlL,'roll'); ...
    getb(mdlL,'steer'); ...
    getb(mdlL,'roll2'); ...
    getb(mdlL,'steer2'); ...
    getb(mdlL,'cross') ];

betaR = [ ...
    getb(mdlR,'(Intercept)'); ...
    getb(mdlR,'roll'); ...
    getb(mdlR,'steer'); ...
    getb(mdlR,'roll2'); ...
    getb(mdlR,'steer2'); ...
    getb(mdlR,'cross') ];

save camber_models_fast.mat betaL betaR
