function car = events_pm(car)
    autox_track = track_pm.loadFromMat('michigantrack2024.mat');
    endurance_track = track_pm.loadFromMat('2024endurancetrack.mat');
    vmax = 30;
    sim = lapsim_pm();
    
    out = struct();
    out.carObj = car.carObj;
    out.spec = car.spec;
    out.comp = struct();
    out.comp.times = struct();
    out.comp.autox = struct();
    out.comp.endurance = struct();
    out.comp.skidpad = struct();
    out.comp.accel = struct();

    %% autox
    [lapT, vprof, vlim] = sim.run(car.carObj, autox_track, vmax);
    out.times.autoX  = lapT;
    out.autox.v_profile = vprof;
    out.autox.v_latlim  = vlim;
    
    %% endurance
    [lapT_endurance, vprof_endurance, vlim_endurance] = sim.run(car.carObj, endurance_track, vmax);
    out.times.endurance = lapT_endurance;
    out.endurance.v_profile = vprof_endurance;
    out.endurance.v_latlim  = vlim_endurance;
    
    %% skidpad
    skidpad = sim.skidpad_pm(car.carObj, 7.625 + car.spec.track_width / 2);
    out.comp.skidpad = skidpad;
    out.times.skidpad = skidpad.time;

    %% accel
    accel = sim.accel(car.carObj);
    out.comp.accel = accel;
    out.times.accel = accel.time;

    car = out;
end