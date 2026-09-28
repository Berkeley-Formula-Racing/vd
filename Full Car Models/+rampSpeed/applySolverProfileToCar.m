function car = applySolverProfileToCar(car,profile)
%APPLYSOLVERPROFILETOCAR Apply profile aero selection to a local Car value.

profile = rampSpeed.resolveSolverProfile(profile);
if profile.aeroMode ~= "static"
    return
end
if ~isobject(car) || ~isprop(car,"rideHeightAero")
    error("rampSpeed:invalidCarAeroConfiguration", ...
        "Approximate aero preview requires a Car with rideHeightAero configuration.");
end

rideHeightAero = car.rideHeightAero;
if ~isstruct(rideHeightAero) || ~isfield(rideHeightAero,"enabled")
    error("rampSpeed:invalidCarAeroConfiguration", ...
        "Approximate aero preview requires rideHeightAero.enabled.");
end
rideHeightAero.enabled = false;
car.rideHeightAero = rideHeightAero;
end
