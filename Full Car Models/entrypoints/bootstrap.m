% BOOTSTRAP Configure paths and the working directory for a Full Car Models entrypoint.
% This script is intended to be run by the launcher scripts in this folder.

entrypointDir = fileparts(mfilename('fullpath'));
fullCarModelsDir = fileparts(entrypointDir);
repoRoot = fileparts(fullCarModelsDir);

addpath(genpath(fullCarModelsDir));
cd(repoRoot);
setup_paths;
