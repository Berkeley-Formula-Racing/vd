function torque_fn = TorqueSpline(tq0, tqL, tqP, tqH)

if nargin < 4
    tq0 = 0;
    tqL = 20.4;
    tqP = 30.5;
    tqH = 10;
end

params.r0           = 0;          
params.tq0          = tq0;          
params.rL           = 5400;       % low→mid boundary RPM
params.tqL          = tqL;       % torque at rL
params.rPeak        = 8500;      % midrange peak RPM
params.tqPeak       = tqP;      % torque at midrange peak
params.rH           = 12000;      % high-end boundary RPM
params.tqH          = tqH;       % torque at rH
params.mid_linearity= 0.65;       % 0 = straight line mid, 1 = very bowed

[rFine_param, tqFine_param, eval_param] = make_torque_curve(params, [0, 12000]);

torque_fn = [rFine_param; tqFine_param];

figure;hold on; grid on; box on
plot(rFine_param, tqFine_param, '-', 'LineWidth',1.8, 'DisplayName','Param curve');
plot([params.r0, params.rL, params.rPeak, params.rH], ...
     [params.tq0, params.tqL, params.tqPeak, params.tqH], 'o', ...
     'DisplayName','Anchors');
xline(params.rL, '--', 'Low→Mid','HandleVisibility','off');
xline(params.rPeak, ':', 'Peak','HandleVisibility','off');
xline(params.rH, '--', 'High boundary','HandleVisibility','off');
xlabel('RPM'); ylabel('Torque');
title('Parametrized Torque Curve (3 segments)');
legend('Location','best');


function [rFine, tqFine, eval_fn] = make_torque_curve(p, rRange)
    assert(p.r0 <= p.rL && p.rL < p.rPeak && p.rPeak <= p.rH, ...
        'Must satisfy r0 ≤ rL < rPeak ≤ rH');

    N = max(2000, round((rRange(2)-rRange(1))/5));
    rFine = linspace(rRange(1), rRange(2), N).';
    anchors = unique([p.r0; p.rL; p.rPeak; p.rH]);
    rFine = unique([rFine; anchors]);  

    isLow  = rFine <= p.rL;
    tqFine = zeros(size(rFine));
    tqFine(isLow) = lerp_line(rFine(isLow), p.r0, p.tq0, p.rL, p.tqL);

    cx = 0.5*(p.rL + p.rPeak);          
    line_mid_tq = 0.5*(p.tqL + p.tqPeak);
    bow = (p.tqPeak - p.tqL) * (0.65 * p.mid_linearity);
    cy  = line_mid_tq + bow;
   

    isBez  = rFine >  p.rL & rFine <= p.rH;
    x0 = p.rL;  y0 = p.tqL;
    x3 = p.rH;  y3 = p.tqH;          
    dx = x3 - x0;
    t  = (rFine(isBez) - x0) / max(dx, eps);
    t  = max(0, min(1, t));         
    
    t_p = (p.rPeak - x0) / max(dx, eps);
    t_p = max(0, min(1, t_p));
    
    alpha = max(0, min(1, p.mid_linearity));  


    x1 = x0 + dx/2;

    y1c = y0 + 0.5*(y3 - y0);

    b0 = (1-t_p)^2; b1 = 2*(1-t_p)*t_p; b2 = t_p^2;

    den = b1;
    if abs(den) < 1e-12, den = sign(den)*1e-12 + (den==0)*1e-12; end
    s   = (p.tqPeak - (b0*y0 + b1*y1c + b2*y3)) / den;
    y1  = y1c + s;

    B0 = (1-t).^2; B1 = 2*(1-t).*t; B2 = t.^2;
    tqFine(isBez) = B0*y0 + B1*y1 + B2*y3;



    isHigh = rFine > p.rH;
    if any(isHigh)
        tqFine(isHigh) = p.tqH;
    end
    
    tqFine(rFine==p.r0)    = p.tq0;
    tqFine(rFine==p.rL)    = p.tqL;
    tqFine(rFine==p.rPeak) = p.tqPeak;
    tqFine(rFine==p.rH)    = p.tqH;
        eval_fn = @(rpm) eval_param_curve(rpm, p, cx, cy);
end
    
function y = lerp_line(x, x0, y0, x1, y1)
    t = (x - x0) / max(x1 - x0, eps);
    y = y0 + t .* (y1 - y0);
end
    
function tq = eval_param_curve(rpm, p)
    r = rpm(:);
    tq = zeros(size(r));

    isLow = r <= p.rL;
    tq(isLow) = lerp_line(r(isLow), p.r0, p.tq0, p.rL, p.tqL);

    isBez = r > p.rL & r <= p.rH;
    if any(isBez)
        x0 = p.rL;  x3 = p.rH;  dx = x3 - x0;
        y0 = p.tqL; y3 = p.tqH;

        t_p = (p.rPeak - x0)/max(dx,eps);
        b0 = (1-t_p)^3; b1 = 3*(1-t_p)^2*t_p; b2 = 3*(1-t_p)*t_p^2; b3 = t_p^3;
        y1c = y0 + (y3 - y0)/3; y2c = y0 + 2*(y3 - y0)/3;
        lambda = 2*p.mid_linearity; if abs(3*(1-t_p)^2*t_p - lambda*3*(1-t_p)*t_p^2) < 1e-9, lambda = 1; end
        den = (b1 - lambda*b2); s = (p.tqPeak - (b0*y0 + b1*y1c + b2*y2c + b3*y3)) / den;
        y1 = y1c + s; y2 = y2c - lambda*s;

        t = (r(isBez) - x0) / max(dx,eps); t = max(0,min(1,t));
        B0 = (1-t).^3; B1 = 3*(1-t).^2 .* t; B2 = 3*(1-t) .* t.^2; B3 = t.^3;
        tq(isBez) = B0*y0 + B1*y1 + B2*y2 + B3*y3;
    end

    isHigh = r > p.rH;
    if any(isHigh)
        tq(isHigh) = lerp_line(r(isHigh), p.rH, p.tqH, p.rH + 2000, p.tqH);
    end
end

end
