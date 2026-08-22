function tbl = ggMetrics(car,caseLabel)
% Builds a named metric table with one row per solved g-g operating point.
%
% Replaces the positional 31-column x_ss/x_accel/x_braking vectors (see
% generate_table.m) with something you can index by name and filter on.
%
% inputs  car:       a car that has already been through makeGG
%         caseLabel: optional string tagging the row (e.g. a sweep case name)
% output  tbl:       table, one row per operating point. Key columns:
%                    point_type  maxLat | accel | brake
%                    vCar, gLat, gLong, gMag, theta   <- sensitivity axes
%                    exitflag                          <- filter on this
%                    plus loads, slips, cambers, attitude and aero terms
%
% theta is the g-g polar angle in degrees: 0 = pure acceleration,
% +90 = pure cornering, 180 = pure braking.
%
% Note: only the solved (+lateral) points are recorded. makeGG mirrors each
% point to -lat by symmetry; those are not re-solved so they carry no extra
% information and would double-count in any distribution or fit.

if nargin < 2
    caseLabel = "";
end

if isempty(car.ss_info)
    error('ggMetrics:noResults', ...
        'car has no g-g results -- run gg2 and makeGG first');
end

sets = { ...
    'maxLat', dedupeRows(car.ss_info); ...
    'accel',  car.accel_info; ...
    'brake',  car.decel_info};

recs = [];
for s = 1:size(sets,1)
    type = sets{s,1};
    info = sets{s,2};
    for i = 1:size(info,1)
        P = info(i,4:12);          % columns 4:12 are the 9-element state vector
        m = car.metrics(P);
        m.point_type = string(type);
        m.case_label = string(caseLabel);
        m.exitflag   = info(i,1);
        recs = [recs; m]; %#ok<AGROW>
    end
end

tbl = struct2table(recs);

% put the identifying and axis columns first
lead = {'case_label','point_type','exitflag','vCar','gLat','gLong','gMag','theta'};
tbl = movevars(tbl,lead,'Before',1);
end

function out = dedupeRows(info)
% ss_info carries one identical maxLat row per lateral-accel grid column
% (makeGG writes it once per ParamSet), so collapse to unique rows
if isempty(info)
    out = info;
else
    out = unique(info,'rows','stable');
end
end
