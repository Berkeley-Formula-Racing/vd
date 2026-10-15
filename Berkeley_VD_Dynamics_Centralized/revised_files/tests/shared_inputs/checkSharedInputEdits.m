function checkSharedInputEdits(root)
% Edit temporary copies to verify both builders' input sources follow changes.
originalPath = path;
tmp = tempname;
configDir = fullfile(tmp,'Full Car Models','config');
componentDir = fullfile(tmp,'Full Car Models','carComponents');
mkdir(configDir); mkdir(componentDir);
cleanup = onCleanup(@() restore(originalPath,tmp)); %#ok<NASGU>
copyfile(fullfile(root,'config','*.m'),configDir);
t = vdTireConfig();
for name = {'Fx_file','Fy_file','camber_ratio_file','camber_model_file'}
    copyfile(t.(name{1}),componentDir);
end
rewrite(fullfile(configDir,'vdVehicleConfig.m'),'carParams.mass = 162;', ...
    'carParams.mass = 172;');
rewrite(fullfile(configDir,'vdTireConfig.m'),'cfg.p_i = 12;', ...
    'cfg.p_i = 13;');
addpath(configDir,'-begin');
clear vdVehicleConfig vdTireConfig vdRampConfig
v = vdVehicleConfig(); t = vdTireConfig(); r = vdRampConfig();
assert(v.car.mass==172 && r.fixedParameters.vehicle.mass==172);
assert(t.p_i==13 && r.fixedParameters.vehicle.p_i==13);
assert(t.studies.camber2_Testcamber.p_i==13);
% New reads reflect another edit, with no persistent configuration cache.
rewrite(fullfile(configDir,'vdVehicleConfig.m'),'carParams.mass = 172;', ...
    'carParams.mass = 173;');
clear vdVehicleConfig vdRampConfig
v = vdVehicleConfig(); r = vdRampConfig();
assert(v.car.mass==173 && r.fixedParameters.vehicle.mass==173);
% Path/settings consumers must not need unrelated full-car coefficient files.
delete(t.Fx_file); delete(t.Fy_file);
pathsOnly = vdTireConfig([],[],false);
assert(pathsOnly.p_i==13 && ~isfield(pathsOnly,'Fx_parameters'));
end
function rewrite(file,old,new)
s = fileread(file); assert(~isempty(strfind(s,old))); %#ok<STREMP>
f = fopen(file,'w'); assert(f>=0);
cleanup = onCleanup(@() fclose(f)); %#ok<NASGU>
fprintf(f,'%s',strrep(s,old,new));
end
function restore(originalPath,tmp)
path(originalPath);
clear vdVehicleConfig vdTireConfig vdRampConfig
if exist(tmp,'dir'), rmdir(tmp,'s'); end
end
