function snapshot = carConfigurationSnapshot(value)
%CARCONFIGURATIONSNAPSHOT Preserve public constructor fields for a car graph.
% Numeric, string, table, and other native data retain their original classes.
snapshot = captureValue(value);
end

function output = captureValue(value)
if isa(value,'scatteredInterpolant')
    output = struct('class',class(value), ...
        'propertyNames',{{'Points';'Values';'Method';'ExtrapolationMethod'}}, ...
        'values',struct('Points',value.Points,'Values',value.Values, ...
        'Method',value.Method,'ExtrapolationMethod',value.ExtrapolationMethod));
elseif isa(value,'griddedInterpolant')
    output = struct('class',class(value), ...
        'propertyNames',{{'GridVectors';'Values';'Method';'ExtrapolationMethod'}}, ...
        'values',struct('GridVectors',{value.GridVectors},'Values',value.Values, ...
        'Method',value.Method,'ExtrapolationMethod',value.ExtrapolationMethod));
elseif isa(value,'triangulation')
    output = struct('class',class(value), ...
        'propertyNames',{{'Points';'ConnectivityList'}}, ...
        'values',struct('Points',value.Points, ...
        'ConnectivityList',value.ConnectivityList));
elseif isa(value,'function_handle')
    output = struct('functionText',func2str(value));
elseif isstruct(value)
    output = value;
    fields = fieldnames(value);
    for element = 1:numel(value)
        for k = 1:numel(fields)
            name = fields{k};
            output(element).(name) = captureValue(value(element).(name));
        end
    end
elseif iscell(value)
    output = value;
    for k = 1:numel(value)
        output{k} = captureValue(value{k});
    end
elseif isobject(value) && any(strcmp(class(value), ...
        {'Car','Aero','Powertrain','Tire2','NondimensionalTire','AeroMap'}))
    names = properties(value);
    output = struct('class',class(value), ...
        'propertyNames',{names},'values',struct());
    if ~isscalar(value)
        output.values = cell(size(value));
        for element = 1:numel(value)
            output.values{element} = captureObject(value(element),names);
        end
    else
        output.values = captureObject(value,names);
    end
else
    output = value;
end
end

function values = captureObject(owner,names)
values = struct();
for k = 1:numel(names)
    name = names{k};
    values.(name) = captureValue(owner.(name));
end
end
