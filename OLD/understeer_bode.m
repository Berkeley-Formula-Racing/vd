%% Bode plots from understeer gradient (deg/g @ tire)
% Requires Control System Toolbox.

clear; clc;

%% ------------ User inputs ------------
% Geometry / inertia / mass
L   =  62 * 0.0254;        % wheelbase [m]
a   = L * 0.512;        % CG -> front axle [m] 
m   = 535 * 0.453592;   % mass [kg] 
Iz  = 83.28;       % yaw inertia [kg*m^2] 

% Understeer gradient at the tire
K_deg_per_g = 0.18;    % [deg/g]  <-- set to your value

% Cornering “size” (only needed to pin down Cf, Cr distribution)
Csum = 1.6e5;      % Cf + Cr [N/rad] (pick based on tire data; adjust as needed)

% Speed sweep for Bode overlays
U_vec = [10 15 20];    % [m/s] add more speeds if you like

% Frequency vector (optional, to standardize all plots)
w = logspace(-1, 2, 400);  % [rad/s]
%% ------------------------------------

g  = 9.80665;
b  = L - a;
Ku = (K_deg_per_g * pi/180) / g;   % SI understeer coeff: rad / (m/s^2)

% Pre-create figure axes
f1 = figure('Name','Yaw-rate / Steer Bode');  ax1 = gca;  hold(ax1,'on');
f2 = figure('Name','Lat-Accel / Steer Bode'); ax2 = gca;  hold(ax2,'on');

leg_txt = strings(0);

for Ui = U_vec(:).'
    % ----- Solve for Cf, Cr that satisfy Ku & Csum -----
    % Equations: Cf + Cr = Csum
    %            (a/Cf - b/Cr) = Ku / m
    % -> quadratic in Cf:  (-Ku)*x^2 + (Ku*Csum + m*L)*x - m*a*Csum = 0
    A = -Ku;
    B = Ku*Csum + m*L;
    C = -m*a*Csum;
    disc = B.^2 - 4*A*C;

    if disc < 0
        error('No real (Cf,Cr) solution for these inputs. Adjust Csum or Ku.');
    end

    Cf_candidates = [(-B + sqrt(disc))/(2*A), (-B - sqrt(disc))/(2*A)];
    % Choose physical root: 0 < Cf < Csum and Cr = Csum - Cf > 0
    Cf = NaN; Cr = NaN;
    for c = Cf_candidates
        if c > 0 && c < Csum
            cr_try = Csum - c;
            % Also prefer positive “understeer-ish” load distribution (optional)
            Cf = c; Cr = cr_try; break;
        end
    end
    if isnan(Cf)
        error('Could not find physical Cf/Cr. Try a different Csum.');
    end

    % ----- Continuous-time bicycle model (states [v; r]) -----
    % v = lateral velocity at CG [m/s], r = yaw rate [rad/s]
    % Input: delta = tire steer angle [rad]
    % Dynamics (linearized):
    %   [v_dot] = [ -(Cf+Cr)/(mU)         -(U + (a*Cf - b*Cr)/(mU)) ] [v] + [ Cf/m      ] delta
    %   [r_dot]   [ -(a*Cf - b*Cr)/(Iz*U)   -(a^2*Cf + b^2*Cr)/(Iz*U)] [r]   [ a*Cf/Iz  ]
    U  = Ui;
    A_m = [ -(Cf+Cr)/(m*U),          -(U + (a*Cf - b*Cr)/(m*U));
            -(a*Cf - b*Cr)/(Iz*U),   -(a^2*Cf + b^2*Cr)/(Iz*U) ];
    B_m = [ Cf/m;  a*Cf/Iz ];

    % Outputs:
    % 1) yaw rate r
    C_r = [0 1]; D_r = 0;

    % 2) lateral acceleration a_y
    % Using: a_y ≈ (Cf+Cr)/(m*U) * v + (a*Cf - b*Cr)/(m*U) * r - (Cf/m)*delta  (from linear tire)
    C_ay = [(Cf+Cr)/(m*U), (a*Cf - b*Cr)/(m*U)];
    D_ay = -Cf/m;

    sys_r  = ss(A_m, B_m, C_r,  D_r);
    sys_ay = ss(A_m, B_m, C_ay, D_ay);

    % ----- Bode plots (magnitude & phase) -----
    figure(f1); bodemag(sys_r, w); hold on;  % yaw-rate/steer
    figure(f2); bodemag(sys_ay, w); hold on; % ay/steer

    leg_txt(end+1) = sprintf('U = %.0f m/s', U);
end

% Beautify plots & legends
figure(f1);
grid on; legend(leg_txt, 'Location', 'southwest'); title('Yaw-rate / Steer (r/\delta)');

figure(f2);
grid on; legend(leg_txt, 'Location', 'southwest'); title('Lateral Accel / Steer (a_y/\delta)');

%% -------- Optional: show steady-state gains vs K check --------
% For verification at each U, steady-state r/δ gain should ~ U/(L + Ku U^2)
% You can quickly print or compare against DC gain:
if true
    fprintf('\nSteady-state yaw-rate gain (theory vs dcgain):\n');
    for Ui = U_vec(:).'
        U = Ui;  % recompute Cf/Cr for this U as above to keep code compact:
        A = -Ku; B = Ku*Csum + m*L; Cc = -m*a*Csum;
        disc = B.^2 - 4*A*Cc; Cf = (-B + sqrt(disc))/(2*A); Cr = Csum - Cf;
        A_m = [ -(Cf+Cr)/(m*U),          -(U + (a*Cf - b*Cr)/(m*U));
                -(a*Cf - b*Cr)/(Iz*U),   -(a^2*Cf + b^2*Cr)/(Iz*U) ];
        B_m = [ Cf/m;  a*Cf/Iz ];
        C_r = [0 1]; D_r = 0;
        sys_r  = ss(A_m, B_m, C_r,  D_r);
        Gdc = dcgain(sys_r);
        Gth = U / (L + Ku*U^2);
        fprintf('U=%4.0f m/s: theory=%.4f, dcgain=%.4f\n', U, Gth, Gdc);
    end
end
