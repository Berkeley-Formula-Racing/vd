function displayed = displayUnits(values,fromUnits,toUnits)
%DISPLAYUNITS Convert values for presentation without changing stored SI data.

if nargin ~= 3
    error('rampSpeed:displayUnits:InvalidArguments', ...
        'displayUnits requires values, fromUnits, and toUnits.');
end

from = normalizeUnit(fromUnits);
to = normalizeUnit(toUnits);
if from == to
    displayed = values;
    return
end

factor = conversionFactor(from,to);
displayed = values .* factor;
end

function unit = normalizeUnit(value)
unit = lower(strtrim(string(value)));
if ~isscalar(unit)
    error('rampSpeed:displayUnits:InvalidUnit','Units must be scalar text.');
end

unit = replace(unit,[char(178) char(179)],['2' '3']);
unit = replace(unit,["²","³"],["2","3"]);
unit = replace(unit,[" ","_","-"],"");
switch unit
    case {"m/s","mps","meterpersecond","meterspersecond"}
        unit = "m/s";
    case {"mph","mi/h","milesperhour"}
        unit = "mph";
    case {"g","gravity"}
        unit = "g";
    case {"m/s2","mps2","m/s^2","mps^2","meterpersecond2"}
        unit = "m/s^2";
    case {"n","newton","newtons"}
        unit = "N";
    case {"lbf","lbforce","poundforce"}
        unit = "lbf";
    case {"m","meter","meters"}
        unit = "m";
    case {"in","inch","inches"}
        unit = "in";
    case {"ft","foot","feet"}
        unit = "ft";
    case {"rad","radian","radians"}
        unit = "rad";
    case {"deg","degree","degrees"}
        unit = "deg";
    case {"rad/(m/s2)","radpermps2","radperm/s2"}
        unit = "rad/(m/s^2)";
    case {"deg/g","degg","degperg"}
        unit = "deg/g";
    case {"rad/s","rps","radpers"}
        unit = "rad/s";
    otherwise
        error('rampSpeed:displayUnits:UnsupportedUnit', ...
            'Unsupported unit "%s".',unit);
end
end

function factor = conversionFactor(from,to)
g = 9.80665;
switch from + "->" + to
    case "m/s->mph"
        factor = 2.2369362920544;
    case "mph->m/s"
        factor = 0.44704;
    case "g->m/s^2"
        factor = g;
    case "m/s^2->g"
        factor = 1/g;
    case "N->lbf"
        factor = 1/4.4482216152605;
    case "lbf->N"
        factor = 4.4482216152605;
    case "m->in"
        factor = 39.37007874015748;
    case "in->m"
        factor = 0.0254;
    case "m->ft"
        factor = 3.280839895013123;
    case "ft->m"
        factor = 0.3048;
    case "rad->deg"
        factor = 180/pi;
    case "deg->rad"
        factor = pi/180;
    case "rad/(m/s^2)->deg/g"
        factor = g*180/pi;
    case "deg/g->rad/(m/s^2)"
        factor = pi/180/g;
    case "rad/s->rps"
        factor = 1/(2*pi);
    case "rps->rad/s"
        factor = 2*pi;
    otherwise
        error('rampSpeed:displayUnits:IncompatibleUnits', ...
            'No conversion is defined from "%s" to "%s".',from,to);
end
end
