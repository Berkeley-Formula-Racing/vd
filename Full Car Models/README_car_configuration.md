# Car configuration and rebuilding

The supported editable source for a QSS car is the scalar setup returned by
`vehicleConfig.defaultSetup`. Change that setup, then call
`vehicleConfig.buildCar` to construct a fresh `Car`. The builder recalculates
derived fields from the edited values and does not carry over a prior solve or
event result.

```matlab
setup = vehicleConfig.defaultSetup("legacy","nominal");
setup.id = "autocross-light";
setup.label = "Autocross light";
setup.parameters.mass = 160;                    % kg
setup.parameters.spring_rate_front_lb_in = 325; % lb/in

[lapCar,lapDerived] = vehicleConfig.buildCar(setup,"lap");
[accelCar,accelDerived] = vehicleConfig.buildCar(setup,"accel");
```

The setup is a value struct. Assigning it and changing one field leaves the
original setup unchanged. `Car`, `Aero`, `Powertrain`, `Tire2`, and
`NondimensionalTire` are value classes too, so assigning a car does not create
a shared mutable setup. A rebuilt car is a new configuration and starts without
the previous car's G-G data or event-result objects. Rebuild after editing;
changing a setup variable does not update a car that was already constructed.

`carConfig` keeps its legacy cell-array result by default. It returns an
N-by-2 cell array, with lap cars in column 1 and acceleration cars in column 2.
The zero-to-four-input calls and one-to-four-output forms remain available.
Use the fifth `"cases"` input to receive named cases with their source setup and
both cars:

```matlab
cars = carConfig();
lapCar = cars{1,1};
accelCar = cars{1,2};

namedCases = carConfig("FullFactorial",[],"legacy","nominal","cases");
firstLapCar = namedCases(1).lapCar;
```

For an explicit design, provide a table whose variable names are fields from
`vehicleConfig.parameterNames`:

```matlab
overrides = table([160;165],"VariableNames",{"mass"});
namedCases = carConfig("Explicit",overrides, ...
    "legacy","nominal","cases");
```

The default values are assembled in `vehicleConfig.defaultSetup`; they are not
stored in `carConfig.m`. A one-off change belongs in the setup or Explicit
design passed to the builder. `README_ramp_speed_app.md` describes the separate
editable setup workflow used by the Ramp Speed app.

## Units and map behavior

QSS setup fields retain the units used by the existing vehicle models. Mass and
weights are in kg; vehicle lengths such as wheelbase, track width, wheel radius,
and CG height are in m. Spring rates use lb/in and motion ratios are
unitless. Tire pressure `p_i` uses psi. Ride heights and the aero map's front
and rear reference heights use inches. Camber angles and the named
`*_deg_per_*` rates use degrees. Suspension ARB stiffness is in Nm/rad; derived
wheel rates are in N/m. Ramp result tables and exports use canonical SI units
such as m/s, m/s^2, N, m, rad, and rad/s.

The QSS setup has its own aero-map path, enable flag, coefficient scales,
centre-of-pressure correction, and reference ride heights. When mapping is
enabled, QSS lap-car construction uses the selected map; acceleration-car
construction uses the configured static acceleration aero values.

Ramp Speed has an independent frozen baseline created by
`carConfigBaseline()`. It does not inherit the QSS design defaults in
`vehicleConfig.defaultSetup`. The current ramp baseline selects `aeromap_b26.csv`
and references its map at the axle-plane extrapolation of the CAD ride-height
line: 3.865589 in front and 5.769078 in rear. The baseline ride heights use the
same datums, giving zero initial map offset. App setup edits rebuild from the
frozen ramp configuration and selected map ID. The app's setup fields use
spring rates in lb/in, ride heights in inches, driver weight in kg, rear weight
distribution in percent, and ARB stiffness in Nm/rad.

## MAT files and saved cars

Existing `carCell` MAT files and ramp study/checkpoint schemas remain supported;
this refactor does not convert those files in place. A saved `Car` is a MATLAB
class object tied to the class definition available when it is loaded, so treat
it as saved simulation state rather than as the editable source of truth. Keep
the source setup or design table when you need to rebuild with current defaults
or current class code. Ramp Speed studies with setup specifications rebuild
cars from those specifications. Older studies without setup specifications
remain readable as result data and are read-only in the app.