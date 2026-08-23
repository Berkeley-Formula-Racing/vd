function [row,details] = doeCaseMetrics(car,caseIndex,rampResult)
%DOECASEMETRICS Reduce one solved DOE car to a design-analysis row.

if nargin < 3, rampResult = []; end
r = blankRow();
r.case_index = caseIndex;
details = struct('ramp',rampResult);

if isempty(car) || isempty(car.comp) || ~isfield(car.comp.times,'autocross')
    row = struct2table(r);
    return
end

c = car.comp;
r.t_autox = c.times.autocross;
r.t_accel = c.times.accel;
r.t_skid = c.times.skidpad;

v = c.autocross.long_vel(:);
ax = c.autocross.long_accel(:);
ay = c.autocross.lat_accel(:);
t = c.autocross.time_vec(:);
dt = [t(1);diff(t)];
good = isfinite(v) & isfinite(ax) & isfinite(ay) & isfinite(dt) & dt >= 0;
power = max(0,car.M*ax(good)).*v(good);
r.total_work_kJ = sum(power.*dt(good))/1000;
r.v_mean_mps = sum(v(good).*dt(good))/max(sum(dt(good)),eps);
r.v_max_mps = max(v(good));
r.gLat_peak_g = max(abs(ay(good)))/car.g;
r.gLat_rms_g = sqrt(sum((ay(good)/car.g).^2.*dt(good))/max(sum(dt(good)),eps));
r.gLong_accel_peak_g = max(ax(good))/car.g;
r.gLong_brake_peak_g = abs(min(ax(good)))/car.g;
r.gLong_rms_g = sqrt(sum((ax(good)/car.g).^2.*dt(good))/max(sum(dt(good)),eps));

G = ggMetrics(car,"case"+caseIndex);
lat = G(G.point_type=="maxLat",:);
acc = G(G.point_type=="accel",:);
brk = G(G.point_type=="brake",:);

speeds = [10 20 30];
for k = 1:numel(speeds)
    tag = string(speeds(k));
    r.("gg_lat_"+tag+"_g") = interpSafe(lat.vCar,abs(lat.gLat),speeds(k));
    r.("gg_accel_"+tag+"_g") = branchAtSpeed(acc,speeds(k),false);
    r.("gg_brake_"+tag+"_g") = branchAtSpeed(brk,speeds(k),true);
end

r.min_Fz_N = min(G.min_Fz,[],'omitnan');
r.wheel_lift_fraction = mean(G.min_Fz <= 0,'omitnan');
r.gg_coverage = coverage(lat,acc,brk);
r.mechanical_balance = median(lat.LLTD,'omitnan');
r.aero_balance = median(lat.CoP,'omitnan');

frontAlpha = mean(abs([lat.alpha_1 lat.alpha_2]),2,'omitnan');
rearAlpha = mean(abs([lat.alpha_3 lat.alpha_4]),2,'omitnan');
proxy = frontAlpha-rearAlpha;
r.understeer_proxy_10_deg = interpSafe(lat.vCar,proxy,10);
r.understeer_proxy_25_deg = interpSafe(lat.vCar,proxy,25);
r.max_constraint_residual = max(abs([lat.lat_accel_residual; ...
    lat.yaw_accel_residual; acc.lat_accel_residual; ...
    acc.yaw_accel_residual; brk.lat_accel_residual; ...
    brk.yaw_accel_residual]),[],'omitnan');

r.valid = all(isfinite([r.t_autox r.t_accel r.t_skid r.total_work_kJ])) ...
    && r.gg_coverage > 0;
if isfinite(r.gg_coverage)
    r.solve_failure_fraction = max(1-r.gg_coverage,0);
end

if ~isempty(rampResult)
    r.understeer_gradient_10_deg_per_g = nearestAtSpeed( ...
        rampResult.perSpeed.vCar,rampResult.perSpeed.K_linear,10);
    r.understeer_gradient_25_deg_per_g = nearestAtSpeed( ...
        rampResult.perSpeed.vCar,rampResult.perSpeed.K_linear,25);
    r.rebalance_speed_mps = firstSignCrossing( ...
        rampResult.perSpeed.vCar,rampResult.perSpeed.K_linear);
end

row = struct2table(r);
end

function y = branchAtSpeed(T,target,asMagnitude)
v = unique(T.vCar);
yv = NaN(size(v));
for i = 1:numel(v)
    rows = find(abs(T.vCar-v(i)) < 1e-9);
    [~,j] = min(abs(T.gLat(rows)));
    yv(i) = T.gLong(rows(j));
end
y = interpSafe(v,yv,target);
if asMagnitude, y = abs(y); end
end

function out = coverage(lat,acc,brk)
v = unique(lat.vCar);
if isempty(v), out = 0; return, end
counts = zeros(numel(v),2);
for i = 1:numel(v)
    counts(i,1) = sum(abs(acc.vCar-v(i)) < 1e-9);
    counts(i,2) = sum(abs(brk.vCar-v(i)) < 1e-9);
end
nLat = max(counts,[],'all');
expected = 2*numel(v)*nLat;
out = min((height(acc)+height(brk))/max(expected,1),1);
end

function y = interpSafe(x,v,xq)
good = isfinite(x) & isfinite(v);
x = x(good); v = v(good);
[x,ia] = unique(x,'stable'); v = v(ia);
[x,order] = sort(x); v = v(order);
if numel(x) < 2 || xq < x(1) || xq > x(end)
    y = NaN;
else
    y = interp1(x,v,xq,'linear');
end
end

function y = nearestAtSpeed(x,v,target)
good = isfinite(x) & isfinite(v);
x = x(good); v = v(good);
if isempty(x)
    y = NaN;
    return
end
[~,i] = min(abs(x-target));
y = v(i);
end

function x0 = firstSignCrossing(x,y)
good = isfinite(x) & isfinite(y);
x = x(good); y = y(good);
[x,order] = sort(x); y = y(order);
for i = 1:numel(x)
    if y(i) == 0
        x0 = x(i);
        return
    end
    if i < numel(x) && sign(y(i)) ~= sign(y(i+1))
        x0 = x(i)-y(i)*(x(i+1)-x(i))/(y(i+1)-y(i));
        return
    end
end
x0 = NaN;
end

function r = blankRow()
r = struct( ...
    'case_index',NaN,'valid',false,'error_message',"", ...
    't_autox',NaN,'t_accel',NaN,'t_skid',NaN,'total_work_kJ',NaN, ...
    'v_mean_mps',NaN,'v_max_mps',NaN,'gLat_peak_g',NaN,'gLat_rms_g',NaN, ...
    'gLong_accel_peak_g',NaN,'gLong_brake_peak_g',NaN,'gLong_rms_g',NaN, ...
    'gg_lat_10_g',NaN,'gg_lat_20_g',NaN,'gg_lat_30_g',NaN, ...
    'gg_accel_10_g',NaN,'gg_accel_20_g',NaN,'gg_accel_30_g',NaN, ...
    'gg_brake_10_g',NaN,'gg_brake_20_g',NaN,'gg_brake_30_g',NaN, ...
    'min_Fz_N',NaN,'wheel_lift_fraction',NaN,'gg_coverage',NaN, ...
    'mechanical_balance',NaN,'aero_balance',NaN, ...
    'understeer_proxy_10_deg',NaN,'understeer_proxy_25_deg',NaN, ...
    'max_constraint_residual',NaN,'solve_failure_fraction',NaN, ...
    'understeer_gradient_10_deg_per_g',NaN, ...
    'understeer_gradient_25_deg_per_g',NaN,'rebalance_speed_mps',NaN);
end
