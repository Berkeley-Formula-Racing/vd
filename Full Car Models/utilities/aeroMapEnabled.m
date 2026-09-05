function tf = aeroMapEnabled(aeroParams)
%AEROMAPENABLED Return whether lap aero should use the ride-height aeromap.
%   use_aeromap is the explicit user-facing switch. map_enabled is retained
%   solely for compatibility with older saved parameter structures.

if isfield(aeroParams,'use_aeromap')
    tf = aeroParams.use_aeromap;
elseif isfield(aeroParams,'map_enabled')
    tf = aeroParams.map_enabled;
else
    tf = false;
end

validateattributes(tf,{'numeric','logical'},{'scalar'},mfilename,'use_aeromap');
tf = logical(tf);
end
