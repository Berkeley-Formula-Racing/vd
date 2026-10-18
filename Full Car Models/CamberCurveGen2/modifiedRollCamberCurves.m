function[camber_matrix] = modifiedRollCamberCurves(roll_angle_vector)
    % Returns phyical camber givn a vector of roll angles
    % Convention: 1 = front left, 2 = front right, 
    %             3 = rear left, 4 = rear right

    % The slopes for linear approximations:
    % taken from WinGeo B20 Final Front Geometry
    ma = 0;
    mb = 0;
    b1 = 1; % static camber, 0 in WinGeo, -1(?) on actual car
    % from B20 Final Rear Geometry
    mc = 0;
    md = 0;
    b2 = 1;

    camber_matrix = zeros(4, numel(roll_angle_vector)); % rows = tire
    for i = 1:numel(roll_angle_vector)
        if roll_angle_vector(i) >= 0
            camber_matrix(1,i) = -mb*roll_angle_vector(i) + b1;
            camber_matrix(2,i) = ma*roll_angle_vector(i) + b1;

            camber_matrix(3,i) = -md*roll_angle_vector(i) + b2;
            camber_matrix(4,i) = mc*roll_angle_vector(i) + b2;

        else if roll_angle_vector(i) < 0
            camber_matrix(1,i) = -ma*roll_angle_vector(i) + b1;
            camber_matrix(2,i) = mb*roll_angle_vector(i) + b1;

            camber_matrix(3,i) = -mc*roll_angle_vector(i) + b2;
            camber_matrix(4,i) = md*roll_angle_vector(i) + b2;
            end
        end

    end
end
