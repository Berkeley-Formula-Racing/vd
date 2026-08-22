classdef Tire2
    properties
        gamma
        p_i
        Fx_parameters
        Fy_parameters
        friction_scaling_factor

        camber2indices
        camber4indices
        camber2ratio
        camber4ratio

    end

    methods
        function obj = Tire2(p_i,Fx_parameters,Fy_parameters,friction_scaling_factor)
            %obj.gamma = gamma;
            obj.p_i = p_i;
            obj.Fx_parameters = Fx_parameters;
            obj.Fy_parameters = Fy_parameters;
            obj.friction_scaling_factor = friction_scaling_factor;
            % named explicitly: a bare load() drops whatever the file holds
            % into this scope and would shadow any same-named function
            c = load("camberratiossmoothed.mat", ...
                'camber2indices','camber4indices','camber2ratio','camber4ratio');
            obj.camber2indices = c.camber2indices;
            obj.camber4indices = c.camber4indices;
            obj.camber2ratio = c.camber2ratio;
            obj.camber4ratio = c.camber4ratio;
        end

        function out = F_y(obj,alpha,kappa,F_z,gamma)
            cambershiftMod = 16.125*gamma*(F_z/250);

            % Camber multiplier: bilinear over (alpha, gamma) on the smoothed
            % ratio tables. This used to be five interp1 calls -- four in
            % alpha, one across gamma -- which cost ~90 us per tire and made
            % F_y 86% of the entire model evaluation. interp1 reparses its
            % inputs and rebuilds an interpolant on every call; a prebuilt
            % griddedInterpolant does the same arithmetic for ~2 us.
            %
            % The 2-D table is mathematically identical to the old chain. The
            % old code interpolated linearly in alpha on each of the four
            % curves, then linearly across gamma -- that IS bilinear, provided
            % the curves are first resampled onto a common alpha axis, and
            % resampling a piecewise-linear function onto a superset of its
            % own breakpoints is exact.
            %
            % Held in a persistent, not a property: the tables come from a
            % fixed file and are the same for every tire, and a new property
            % would break every Tire2 already saved inside a .mat sweep cache.
            persistent camberF warnedNaN
            if isempty(camberF)
                camberF = buildCamberInterp(obj);
            end
            cambermultiplier = camberF(alpha,gamma);

            gamma = 0;%obj.gamma*0.0174533; %degrees to radians

            alpha = alpha*0.0174533; %degrees to radians
            F_z = F_z*0.224809; %N to lbf

            F_z0 = 200; %nominal load
            p_i0 = 13;  %nominal pressure
            F_z = abs(F_z);

            df_z = (F_z- F_z0)./F_z;
            dp_i = (obj.p_i-p_i0)./obj.p_i;

            %% Parameters

            %[p_cy1,p_dy1,p_dy2,p_dy3,p_ey1,p_ey2,p_ey3,p_ey4,p_ey5,p_hy1,p_hy2,p_ky1,p_ky2,p_ky3,...
               % p_ky4,p_ky5,p_ky6,p_ky7,p_py1,p_py2,p_py3,p_py4,p_py5,p_vy1,p_vy2,p_vy3,p_vy4,r_by1,r_by2,...
                %r_by3,r_by4,r_cy1,r_ey1,r_ey2,r_hy1,r_hy2,r_vy1,r_vy2,r_vy3,r_vy4,r_vy5,r_vy6] = obj.Fy_parameters{:};

            % one property read instead of 42: property access on a classdef
            % object is not free, and this runs inside every fmincon iteration
            p = obj.Fy_parameters;

            p_cy1 = p(1);
            p_dy1 = p(2);
            p_dy2 = p(3);
            p_dy3 = p(4);
            p_ey1 = p(5);
            p_ey2 = p(6);
            p_ey3 = p(7);
            p_ey4 = p(8);
            p_ey5 = p(9);
            p_hy1 = p(10);
            p_hy2 = p(11);
            p_ky1 = p(12);
            p_ky2 = p(13);
            p_ky3 = p(14);
            p_ky4 = p(15);
            p_ky5 = p(16);
            p_ky6 = p(17);
            p_ky7 = p(18);
            p_py1 = p(19);
            p_py2 = p(20);
            p_py3 = p(21);
            p_py4 = p(22);
            p_py5 = p(23);
            p_vy1 = p(24);
            p_vy2 = p(25);
            p_vy3 = p(26);
            p_vy4 = p(27);
            r_by1 = p(28);
            r_by2 = p(29);
            r_by3 = p(30);
            r_by4 = p(31);
            r_cy1 = p(32);
            r_ey1 = p(33);
            r_ey2 = p(34);
            r_hy1 = p(35);
            r_hy2 = p(36);
            r_vy1 = p(37);
            r_vy2 = p(38);
            r_vy3 = p(39);
            r_vy4 = p(40);
            r_vy5 = p(41);
            r_vy6 = p(42);

            lambda_cy = 1;       % shape factor
            lambda_ey = 1;       % curvature
            lambda_hy = 0;       % horizontal shift
            lambda_kyalpha = 1;  % cornering stiffness
            lambda_kygamma = 1;  % camber force stiffness
            lambda_muy = 1;      % peak friction coefficient
            lambda_vy = 0;       % vertical shift

            lambda_ykappa = 1;   % influence on F_y(alpha)
            lambda_vykappa = 0;  % induced ply-steer F_y

            %%  Magic Formula Equations
            K_yalpha = p_ky1.*F_z0.*(1+p_py1.*dp_i).*sin(p_ky4.*atan(F_z./...
                ((p_ky2 + p_ky5*gamma.^2).*(1+p_py2*dp_i).*F_z0))).*...
                (1-p_ky3.*abs(gamma)).*lambda_kyalpha;
            K_ygamma = (p_ky6+p_ky7.*df_z).*(1+p_py5.*dp_i).*F_z.*lambda_kygamma;
            S_vy0 = F_z.*(p_vy1+p_vy2.*df_z).*lambda_vy.*lambda_muy;
            S_vygamma = F_z.*(p_vy3+p_vy4.*df_z).*gamma.*lambda_kygamma.*lambda_muy;
            S_vy = S_vy0+S_vygamma;
            S_hy0 = (p_hy1+p_hy2.*df_z).*lambda_hy;
            S_hygamma = (K_ygamma.*gamma-S_vygamma)./K_yalpha;
            S_hy = S_hy0+S_hygamma;
            alpha_y = alpha+S_hy;
            mu_y = (p_dy1+p_dy2.*df_z).*(1-p_dy3.*gamma.^2).*...
                (1+p_py3.*dp_i+p_py4.*dp_i.^2).*lambda_muy;

            C_y = p_cy1.*lambda_cy;
            D_y = mu_y.*F_z;
            E_y = (p_ey1+p_ey2.*df_z).*(1+p_ey5.*gamma.^2-(p_ey3+p_ey4.*gamma).*sign(alpha_y)).*lambda_ey;
            B_y = K_yalpha./C_y./D_y;

            %Combined Slip

            B_ykappa = (r_by1 + r_by4.*gamma.^2).*cos(atan(r_by2.*(alpha - r_by3))).*lambda_ykappa;
            C_ykappa = r_cy1;
            E_ykappa = r_ey1 + r_ey2.*df_z;
            S_hykappa = r_hy1 + r_hy2.*df_z;
            kappa_s = kappa + S_hykappa;
            D_vykappa = mu_y.*F_z.*(r_vy1 + r_vy2.*df_z + r_vy3.*gamma).*cos(atan(r_vy4.*alpha));


            S_vykappa = D_vykappa.*sin(r_vy5.*atan(r_vy6.*kappa)).*lambda_vykappa;
            G_ykappa = cos(C_ykappa.*atan(B_ykappa.*kappa_s - E_ykappa.*(B_ykappa.*kappa_s - ...
                atan(B_ykappa.*kappa_s))))./cos(C_ykappa.*atan(B_ykappa.*S_hykappa - E_ykappa...
                .*(B_ykappa.*S_hykappa - atan(B_ykappa.*S_hykappa))));

            %Lateral Force
            F_y = D_y.*sin(C_y.*atan(B_y.*alpha_y-E_y.*(B_y.*alpha_y-atan(B_y.*alpha_y))))+S_vy;
            F_y = transpose(G_ykappa.*F_y+S_vykappa);

            F_y(F_z==0) = 0; %zero load

            F_y2 = F_y.*4.44822.*obj.friction_scaling_factor; %lbf to N, scaled
            out = F_y2.*cambermultiplier + cambershiftMod;
            % A NaN here means the tire model produced garbage and fmincon is
            % about to search on it. This used to disp() three unlabelled
            % lines with no clue which corner or state caused it -- inside a
            % parfor g-g that is thousands of lines of noise that identify
            % nothing. Report the inputs, once, under an ID that can be muted
            % with warning('off','Tire2:nanForce').
            if isnan(out)
                if isempty(warnedNaN)
                    warnedNaN = true;
                    warning('Tire2:nanForce', ...
                        ['F_y returned NaN: alpha=%g deg, kappa=%g, F_z=%g N, ' ...
                         'gamma=%g deg (pre-camber F_y2=%g). Reported once per ' ...
                         'session.'], alpha/0.0174533, kappa, F_z/0.224809, ...
                         cambershiftMod, F_y2);
                end
            end
        end

        function out = F_x(obj,alpha,kappa,F_z, gamma)

            gamma = gamma*0.0174533; %degrees to radians
            alpha_f = alpha*0.0174533; %degrees to radians
            F_z = F_z*0.224809; %N to lbf

            F_z0 = 200; %nominal load
            p_i0 = 13;  %nominal pressure
            F_z = abs(F_z);

            df_z = (F_z- F_z0)./F_z;
            dp_i = (obj.p_i-p_i0)./obj.p_i;

            %% Parameters

           %[p_cx1,p_dx1,p_dx2,p_dx3,p_ex1,p_ex2,p_ex3,p_ex4,p_hx1,p_hx2,p_kx1,p_kx2,p_kx3,p_px1,...
               % p_px2,p_px3,p_px4,p_vx1,p_vx2,r_bx1,r_bx2,r_bx4,r_cx1,r_ex1,r_ex2,r_hx1] = obj.Fx_parameters{:};

            % one property read instead of 26, as in F_y above
            q = obj.Fx_parameters;

            p_cx1 = q(1);
            p_dx1 = q(2);
            p_dx2 = q(3);
            p_dx3 = q(4);
            p_ex1 = q(5);
            p_ex2 = q(6);
            p_ex3 = q(7);
            p_ex4 = q(8);
            p_hx1 = q(9);
            p_hx2 = q(10);
            p_kx1 = q(11);
            p_kx2 = q(12);
            p_kx3 = q(13);
            p_px1 = q(14);
            p_px2 = q(15);
            p_px3 = q(16);
            p_px4 = q(17);
            p_vx1 = q(18);
            p_vx2 = q(19);
            r_bx1 = q(20);
            r_bx2 = q(21);
            r_bx4 = q(22);
            r_cx1 = q(23);
            r_ex1 = q(24);
            r_ex2 = q(25);
            r_hx1 = q(26);

            lambda_cx  = 1;     %shape factor
            lambda_ex  = 1;     %curve factor
            lambda_hx  = 0;     %horizontal shift
            lambda_kxkappa = 1; %brake slip stiffness
            lambda_mux = 1;      % peak friction coefficient
            lambda_vx  = 0;     %vertical shift

            lambda_xalpha = 1; %influence on F_x(kappa)

            %% Magic Formula Equations
            mu_x = (p_dx1+p_dx2.*df_z).*(1-p_dx3.*gamma.^2).*(1+p_px3.*dp_i+p_px4.*dp_i.^2).*lambda_mux;

            K_xkappa = (p_kx1+p_kx2.*df_z).*exp(p_kx3.*df_z).*(1+p_px1.*dp_i+p_px2.*dp_i.^2).*F_z.*lambda_kxkappa;
            if ~isfinite(K_xkappa)
                K_xkappa = 1e200;
            end
            S_hx = (p_hx1+p_hx2.*df_z).*lambda_hx;
            S_vx = (p_vx1+p_vx2.*df_z).*F_z.*lambda_vx.*lambda_mux;
            kappa_x = kappa + S_hx;
            C_x = p_cx1.*lambda_cx;
            D_x = mu_x.*F_z;
            E_x = (p_ex1+p_ex2.*df_z+p_ex3.*df_z.^2).*(1-p_ex4.*sign(kappa_x)).*lambda_ex;
            B_x = K_xkappa./(C_x.*D_x);

            %Combined slip:
            B_xalpha = (r_bx1+r_bx4.*gamma.^2).*cos(atan(r_bx2.*kappa)).*lambda_xalpha;
            C_xalpha = r_cx1;
            E_xalpha = r_ex1 + r_ex2.*df_z;
            S_hxalpha = r_hx1;
            alpha_s = alpha_f + S_hxalpha;
            G_xalpha = cos(C_xalpha.*atan(B_xalpha.*alpha_s - E_xalpha.*(B_xalpha.*alpha_s - ...
                atan(B_xalpha.*alpha_s))))./cos(C_xalpha.*atan(B_xalpha.*S_hxalpha - E_xalpha...
                .*(B_xalpha.*S_hxalpha - atan(B_xalpha.*S_hxalpha))));

            %Longitudinal Force
            F_x = transpose((D_x.*sin(C_x.*atan(B_x.*kappa_x-E_x.*(B_x.*kappa_x-atan(B_x.*kappa_x))))+S_vx).*G_xalpha);

            F_x(F_z==0) = 0; %zero load

            out = F_x*4.44822*obj.friction_scaling_factor; %lbf to N, scaled

        end

        function F = buildCamberInterp(obj)
            % Collapses the four camber ratio curves into one bilinear
            % interpolant over (alpha, gamma), replacing five interp1 calls
            % per tire evaluation with one lookup.
            %
            % Curve layout, matching what F_y used to compute:
            %   gamma = +4  ->  camber4ratio(alpha)
            %   gamma = +2  ->  camber2ratio(alpha)
            %   gamma =  0  ->  1
            %   gamma = -2  ->  camber2ratio(-alpha)   [the old interp1 over
            %   gamma = -4  ->  camber4ratio(-alpha)    the NEGATED index grid]
            %
            % Every curve is resampled onto the union of all four breakpoint
            % sets. That union contains each curve's own breakpoints, so the
            % resampling is exact for a piecewise-linear function rather than
            % an approximation of it. Linear extrapolation at both ends
            % reproduces interp1's "extrap", including past the short end of
            % the gamma=4 table, since the added node lies on the same line
            % the old code would have extrapolated along.
            a2 = obj.camber2indices(:); r2 = obj.camber2ratio(:);
            a4 = obj.camber4indices(:); r4 = obj.camber4ratio(:);

            aU = unique([a2; a4; -a2; -a4]);

            M = [ interp1(a4,r4,-aU,'linear','extrap'), ...   % gamma = -4
                  interp1(a2,r2,-aU,'linear','extrap'), ...   % gamma = -2
                  ones(numel(aU),1),                    ...   % gamma =  0
                  interp1(a2,r2, aU,'linear','extrap'), ...   % gamma = +2
                  interp1(a4,r4, aU,'linear','extrap') ];     % gamma = +4

            F = griddedInterpolant({aU,[-4 -2 0 2 4]},M,'linear','linear');
        end

    end

end
