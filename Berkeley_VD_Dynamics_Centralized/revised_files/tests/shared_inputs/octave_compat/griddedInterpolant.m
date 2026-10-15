function F = griddedInterpolant(grid,values,method,extrapolation)
% Octave test adapter for the tire's two-dimensional linear interpolation.
% This directory must only be added by runSharedInputChecks in GNU Octave.
assert(iscell(grid) && numel(grid)==2);
assert(strcmp(method,'linear') && strcmp(extrapolation,'linear'));
F = @(x,y) linearGrid(grid{1},grid{2},values,x,y);
end
function out = linearGrid(gx,gy,values,x,y)
if isscalar(x), x = repmat(x,size(y)); end
if isscalar(y), y = repmat(y,size(x)); end
assert(isequal(size(x),size(y)));
v = interp1(gx,values,x(:),'linear','extrap');
pos = interp1(gy,1:numel(gy),y(:),'linear','extrap');
lo = max(1,min(numel(gy)-1,floor(pos)));
w = pos-lo;
i = (1:numel(x))';
out = (1-w).*v(sub2ind(size(v),i,lo)) + w.*v(sub2ind(size(v),i,lo+1));
out = reshape(out,size(x));
end
