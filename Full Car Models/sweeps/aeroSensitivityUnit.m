function unit = aeroSensitivityUnit(S,param,quantity)
%AEROSENSITIVITYUNIT Return a display unit for static or aeromap inputs.

if nargin < 3, quantity = 'g'; end
base = parameterUnit(S,param);
switch lower(quantity)
    case 'time'
        unit = ['s per ' base];
    otherwise
        unit = ['g per ' base];
end
end

function unit = parameterUnit(S,param)
if isfield(S,'paramInfo')
    names = string({S.paramInfo.name});
    index = find(names == string(param),1);
    if ~isempty(index)
        unit = char(S.paramInfo(index).unit);
        return
    end
end
switch char(param)
    case {'ClA','CdA'}, unit = 'm^2';
    case 'CoP', unit = 'unit front balance';
    otherwise, unit = 'unit';
end
end
