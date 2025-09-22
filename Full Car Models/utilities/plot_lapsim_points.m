function plot_lapsim_points(carCell, display_point_values_above_bar_flag, label_cars_automatically_flag,...
    manual_car_labels, automatic_label_name, automatic_label, selected_categories)

    %remove accel cars from carCell
    carCell = carCell(:,1);
    
    % when only one category, use car labels as x labels, when multiple
    % categories, use car label in legend
    single_category_mode = (numel(selected_categories) == 1);

    
    % label cars
    num_cars = numel(carCell);
    if label_cars_automatically_flag
        car_labels = cell(1,num_cars);
        for i = 1:num_cars
            car_labels{i} = char(string(automatic_label(carCell{i})));
        end
    else
        if numel(manual_car_labels) ~= num_cars
            throw(MException('plot_lapsim_points:wrongNumberManualLabels',...
                'number of cars and manual labels must match'));
        end
        car_labels = manual_car_labels;
    end

    points_categories = {'Accel','Autocross','Endurance','Skidpad','Total'};
    point_vector = @(points) [points.accel; points.autocross; points.endurance; points.skidpad; points.total];

    %make array of all cars points in selected categories
    car_point_array = zeros(numel(points_categories),num_cars);
    for i = 1:num_cars
        points = carCell{i}.comp.points;
        point_vector_temp = point_vector(points);
        car_point_array(:,i) = point_vector_temp;
    end
    
    %ensure that matlab does not re-order data
    pre_ordered_categorical = @(x) reordercats(categorical(x), x);

    %plot data
    if single_category_mode
        b = bar(pre_ordered_categorical(car_labels),...
            car_point_array(selected_categories,:));

        title(['Simulated Competition ' points_categories{selected_categories} ' Points']);

        if label_cars_automatically_flag
            xlabel(automatic_label_name);
        end
    else
        b = bar(pre_ordered_categorical(...
            points_categories(selected_categories)),...
            car_point_array(selected_categories,:));

        %legend
        legendObj = legend(car_labels,'location','northwest');
        if label_cars_automatically_flag
            title(legendObj, automatic_label_name);
        end

        title('Simulated Competition Points');
    end

    ylabel('Points');

    %label bar heights
    if display_point_values_above_bar_flag 
        for i = 1:numel(b)
            xtips = b(i).XEndPoints;
            ytips = b(i).YEndPoints;
            labels = string(round(b(i).YData,1));
            text(xtips,ytips,labels,'HorizontalAlignment','center',...
                'VerticalAlignment','bottom','FontSize',7)
        end
    end
end


