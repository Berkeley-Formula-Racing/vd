function grid = ggGrid(maxVelocity,opts)
%GGGRID Select the production or fast-screening g-g resolution.
%   The screening grid is for finding sensitivity trends quickly. Use the
%   production grid for reported performance and lap-time results.

if nargin < 2 || isempty(opts), opts = struct(); end
if ~isstruct(opts) || numel(opts) ~= 1
    error('ggGrid:badOptions','opts must be a scalar struct.');
end
validateattributes(maxVelocity,{'numeric'},{'scalar','real','finite','positive'}, ...
    mfilename,'maxVelocity');

fastScreening = false;
if isfield(opts,'fastScreening')
    fastScreening = opts.fastScreening;
end
validateattributes(fastScreening,{'numeric','logical'}, ...
    {'scalar','real','finite','binary'},mfilename,'opts.fastScreening');

grid.fastScreening = logical(fastScreening);
if grid.fastScreening
    grid.velocityInterval = 2;
    grid.lateralCount = 10;
else
    grid.velocityInterval = 1;
    grid.lateralCount = 20;
end

minVelocity = 5;
grid.velocity = minVelocity:grid.velocityInterval:maxVelocity;
if isempty(grid.velocity)
    grid.velocity = maxVelocity;
elseif maxVelocity-grid.velocity(end) > 0.25*grid.velocityInterval
    grid.velocity(end+1) = maxVelocity;
else
    grid.velocity(end) = maxVelocity;
end
end
