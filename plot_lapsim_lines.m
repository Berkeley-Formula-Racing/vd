function plot_lapsim_lines_by_event(carCell, label_cars_automatically_flag, ...
    manual_car_labels, automatic_label_name, automatic_label)

    % Remove accel cars from carCell
    carCell = carCell(:,1);
    num_cars = numel(carCell);

    % Get labels
    if label_cars_automatically_flag
        car_labels = cell(1,num_cars);
        for i = 1:num_cars
            car_labels{i} = char(string(automatic_label(carCell{i})));
        end
    else
        if numel(manual_car_labels) ~= num_cars
            error('Number of cars and manual labels must match');
        end
        car_labels = manual_car_labels;
    end

    event_names = {'Accel','Autocross','Endurance','Skidpad'};
    num_events = numel(event_names);

    % Build point matrix: rows = events, columns = cars
    points_matrix = zeros(num_events, num_cars);
    for i = 1:num_cars
        pts = carCell{i}.comp.points;
        points_matrix(:,i) = [pts.accel; pts.autocross; pts.endurance; pts.skidpad];
    end

    %% --- Plot: Each line = 1 event, X-axis = cars ---
    figure('Name','Event Points Across Cars'); hold on;
    for e = 1:num_events
        plot(1:num_cars, points_matrix(e,:), '-o', 'DisplayName', event_names{e});
    end
    xticks(1:num_cars);
    xticklabels(car_labels);
    xlabel('Car');
    ylabel('Points');
    title('Points by Event Across Cars');
    legend('Location','northwest');
    grid on;
    hold off;
end
