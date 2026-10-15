function car = vdApplySuspensionConfig(car,study)
%VDAPPLYSUSPENSIONCONFIG Apply a named upstream suspension study setup.
cfg = vdVehicleConfig();
tire = vdTireConfig([],[],false);
if ~isfield(cfg.suspension,study) || strcmp(study,'damping_file')
    error('vdVehicleConfig:unknownStudy','Unknown suspension study: %s',study);
end
s = cfg.suspension.(study);
fields = {'k','Ixx','Iyy','k_rf','k_rr'};
for k = 1:numel(fields), car.(fields{k}) = s.(fields{k}); end
if strcmp(study,'vibration')
    car.c = s.c;
    return
end
if strcmp(study,'transient')
    car.k_tf = tire.vertical_stiffness.transient_N_per_m;
else
    car.k_tf = car.k*tire.vertical_stiffness.modal_spring_multiplier;
end
car.k_tr = car.k_tf;
data = load(cfg.suspension.damping_file,'cdamp','rdamp',s.motion_ratio_front,s.motion_ratio_rear);
car.c_compression = data.cdamp{s.damper_index}([1 2 10:10:end],:);
car.c_rebound = data.rdamp{s.damper_index}([1 2 10:10:end],:);
car.MR_F = data.(s.motion_ratio_front);
car.MR_R = data.(s.motion_ratio_rear);
end
