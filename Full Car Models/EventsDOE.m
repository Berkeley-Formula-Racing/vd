classdef EventsDOE < handle
    % Calculates results for different dynamic events
    
    properties
        car
        accelCar
        
        autocross_track
        
        % info structs
        skidpad
        accel
        autocross
        
        times % Stores the final time for each event in seconds
        
        % contains info used for interpolation in autocross/accel solvers
        interp_info
    end
    
    methods
        function obj = EventsDOE(car,accelCar)
            obj.car = car;
            obj.accelCar = accelCar;
            
            % maps
            load('michigantrack.mat', 'arclength', 'curvature'); 
            obj.autocross_track = [arclength; curvature];
            
            % sweep for max velocity for given radius
            % used for interpolation in Autocross
            [x_table_corner_vel,radius_vector,max_vel_corner_vector] = vel_cornering_sweep(obj.car);
            
            obj.interp_info.x_table_corner_vel = x_table_corner_vel;
            obj.interp_info.radius_vector = radius_vector;
            obj.interp_info.max_vel_corner_vector = max_vel_corner_vector;
            
            % sweep for max pure longitudinal acceleration for given velocity
            % used for interpolation in Accel
            [x_table_accel,long_vel_guess,long_accel_matrix] = long_accel_sweep(obj.accelCar);
            obj.interp_info.x_table_accel = x_table_accel;
            obj.interp_info.long_vel_guess = long_vel_guess;
            obj.interp_info.long_accel_matrix = long_accel_matrix;
        end
        
        function times = calcTimes(obj)
            % Runs all events and returns the raw times struct
            obj.Skidpad();
            obj.Accel();
            obj.Autocross();
            times = obj.times;
        end
        
        function [x_table_skid,maxVel,time] = Skidpad(obj)
            % modelled as pure steady state (no longitudinal acceleration)
            % inner radius of skidpad is 7.625 m
            % width of skidpad is 3 m
            radius = 8.5;
            [x_table_skid,maxVel,time] = max_skidpad_vel(radius,obj.car);
            
            
            obj.times.skidpad = time;
            obj.skidpad.x_table_skid = x_table_skid;
            obj.skidpad.metrics.understeer_deg = x_table_skid.alpha_1 - x_table_skid.alpha_3;
        end
        
        function [time,ending_vel,long_accel_vector,long_vel_vector] = Accel(obj)
            % car starts 0.3 m behind starting line
            % accel is 75 m long
             long_vel = 0;
             long_vel_interp = obj.interp_info.long_vel_guess;
             long_accel_interp = obj.interp_info.long_accel_matrix;
             
             [~,ending_vel,~,~] = straight(long_vel,2.19,long_vel_interp,...
                  long_accel_interp,obj.accelCar.max_vel,obj.accelCar);
              
             % starting velocity for accel is ending velocity of 0.3 straight
             long_vel = ending_vel;
            
            [time_vec,ending_vel,long_accel_vector,long_vel_vector] = straight(long_vel,75,...
                long_vel_interp,long_accel_interp,obj.accelCar.max_vel,obj.accelCar);
            
            obj.times.accel = time_vec(end);
            obj.accel.time_vec = time_vec;
            obj.accel.long_vel_vector = long_vel_vector;
            obj.accel.long_accel_vector = long_accel_vector;
        end
        
        function [long_vel_final,long_accel_final,lat_accel_final,time_final,...
                    time_vec,num_upshifts] = Track_Solver(obj,arclength,curvature, rollout, init_vel)
            % finds apexes in curvature profile (apex of corner) and finds
            % max possible velocity at each apex
            % then max accel and max braking are calculated for each segment
            
            [extrema,extrema_indices] = curvature_apexes(arclength,curvature);
            arclength = [0 arclength];
            
            [apex_velocity] = apex_velocities(obj.interp_info.radius_vector,obj.interp_info.max_vel_corner_vector,extrema);
            
            [F_accel,F_braking] = create_scattered_interpolants2(obj.car.longAccelLookup, obj.car.longDecelLookup);
            
            long_vel = 0;
            if rollout
                long_vel_interp = obj.interp_info.long_vel_guess;
                long_accel_interp = obj.interp_info.long_accel_matrix;
                [~,ending_vel,~,~] = straight(0,6,long_vel_interp,long_accel_interp,obj.car.max_vel,obj.car);
                long_vel = ending_vel;
            else
                long_vel = init_vel;
            end
            
            lat_accel_vector_1 = [];
            long_accel_vector_1 = [];
            long_vel_vector_1 = [];
            time_1 = [];
            lat_accel_vector_2 = [];
            long_accel_vector_2 = [];
            long_vel_vector_2 = [];
            time_2 = [];
            last_index_1 = 1;
            last_index_2 = 1;
            num_upshifts = 0;
            
            for j = 1:numel(extrema_indices) 
                current_gear = find(long_vel<obj.car.powertrain.switch_gear_velocities,1);
                shift_time_cumulative = 0;
                start_shifting = false;
                
                % FORWARD INTEGRATION
                for i = last_index_1:extrema_indices(j)       
                    lat_accel = long_vel^2*abs(curvature(i));
                    lat_accel_vector_1(i) = lat_accel*sign(curvature(i));
                    
                    if long_vel>obj.car.powertrain.switch_gear_velocities(current_gear)
                        start_shifting = true;
                    end
                    
                    if start_shifting && shift_time_cumulative < obj.car.powertrain.shift_time
                        long_accel = -obj.car.aero.drag(long_vel)/obj.car.M;
                        shift_time_cumulative = shift_time_cumulative+time_1(i-1);
                    else
                        long_accel = F_accel(real(lat_accel), real(long_vel));
                    end
                    
                    if shift_time_cumulative>obj.car.powertrain.shift_time
                        start_shifting = false;
                        shift_time_cumulative = 0;
                        current_gear = current_gear+1;
                        num_upshifts = num_upshifts + 1;
                    end
                    
                    long_accel_vector_1(i) = long_accel;
                    long_vel_initial = long_vel;
                    long_vel_vector_1(i) = long_vel_initial;
                    long_vel = sqrt(long_vel^2+2*long_accel*(arclength(i+1)-arclength(i)));
                    time_1(i) = 2*(arclength(i+1)-arclength(i))/(long_vel+long_vel_initial);
                end
                
                last_index_1 = extrema_indices(j);
                long_vel = apex_velocity(j);
                
                % BACKWARD INTEGRATION
                for i = extrema_indices(j):-1:last_index_2 
                    lat_accel = long_vel^2*abs(curvature(i));
                    lat_accel_vector_2(i) = lat_accel*sign(curvature(i));
                    long_accel = F_braking(lat_accel,long_vel);
                    
                    if long_vel == obj.car.max_vel
                        long_accel = 0;
                    end
                    
                    long_accel_vector_2(i) = long_accel;
                    long_vel_initial = long_vel;
                    long_vel_vector_2(i) = long_vel_initial;
                    long_vel = sqrt(long_vel^2-2*long_accel*(arclength(i+1)-arclength(i)));
                    time_2(i) = 2*(arclength(i+1)-arclength(i))/(long_vel+long_vel_initial);
                end
                
                last_index_2 = extrema_indices(j);
                long_vel = min(long_vel_vector_1(end),long_vel_vector_2(end));
            end  
            
            % MINIMUM OF PROFILES
            long_vel_final = min(long_vel_vector_1,long_vel_vector_2);
            indices_2 = find(long_vel_vector_2 == long_vel_final);
            
            long_accel_final = long_accel_vector_1;
            long_accel_final(indices_2) = long_accel_vector_2(indices_2);
            
            lat_accel_final = lat_accel_vector_1;
            lat_accel_final(indices_2) = lat_accel_vector_2(indices_2);
            
            time_final = time_1;
            time_final(indices_2) = time_2(indices_2);
            time_vec = cumsum(time_final);
            time_final = sum(time_final);
        end
        
        function [long_vel_final,long_accel_final,lat_accel_final,time_final] = Autocross(obj)
            arclength = obj.autocross_track(1,:);
            curvature = obj.autocross_track(2,:);
            [long_vel_final,long_accel_final,lat_accel_final,time_final,time_vec,num_upshifts] = ...
                Track_Solver(obj,arclength,curvature, true, 0);
            
            obj.times.autocross = time_final;
            obj.autocross.time_vec = time_vec;
            obj.autocross.num_upshifts = num_upshifts;
            obj.autocross.long_vel = long_vel_final;
            obj.autocross.long_accel = long_accel_final;
            obj.autocross.lat_accel = lat_accel_final;
            obj.autocross.metrics.autox_time = time_final;
            obj.autocross.metrics.total_work_J = sum( max(0, long_accel_final .* obj.car.M) .* long_vel_final .* time_final );
            obj.autocross.metrics.max_lat_g = max(abs(lat_accel_final)) / 9.81;
        end
    end
end