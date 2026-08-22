folder = fileparts(which('setup_paths'));
addpath(genpath(folder));

% Warn when something OUTSIDE this folder is shadowing a file inside it.
%
% addpath prepends, so this folder wins at the moment setup_paths runs -- but
% any later addpath of a sibling directory silently takes precedence, and the
% repo has real duplicates: Other/BDI Tire Model/ holds a Tire2 whose
% constructor takes gamma as an extra first argument. Loading that one instead
% is not a subtle failure to debug from the symptom; parameters_loop passes 4
% arguments to a 5-argument constructor and the error lands on a line number
% inside a file the caller never mentions. Had the signatures matched it would
% not have errored at all -- a different tire model would just have been used.
%
% Checked for the classes and functions the lap sim actually calls, not
% everything: a duplicate somewhere unused is not worth a warning at startup.
checkShadowing(folder);

% Project plot font, set on the graphics root so EVERY figure created from
% here on inherits it -- including the plotting functions that never ask for
% it. Must come after addpath, since setPlotFont lives in sweeps/.
% Warns once per session if the font is not visible to MATLAB.
setPlotFont();




function checkShadowing(folder)
% Two different failures, one message.
%
% 1. PATH SHADOW -- a copy outside this folder wins on the path. addpath
%    prepends, so this folder wins the moment setup_paths runs, but any later
%    addpath of a sibling takes precedence silently.
% 2. STALE CLASS -- MATLAB caches a classdef once it is loaded, so if the wrong
%    copy was loaded earlier in the session the PATH can be correct while the
%    LOADED definition is still the other one. "clear" does not clear classes.
%    This is the one that actually bites, and it presents as a
%    wrong-argument-count error inside a file the caller never mentions.
%
% The repo really does have duplicates: Other/BDI Tire Model/ holds a Tire2
% whose constructor takes gamma as an extra first argument. parameters_loop
% passes 4 arguments to that 5-argument constructor and the error surfaces deep
% inside Tire2. Had the signatures matched it would not have errored at all --
% a different tyre model would simply have been used.
%
% Arity is how case 2 is caught: these are the project's constructors, and a
% mismatch means something else is loaded under that name.
watch = {'Tire2','Car','Aero','Powertrain','Events2','gg2','makeGG', ...
         'carConfig','parameters_loop','straight','max_lat_accel'};
arity = struct('Tire2',4,'Aero',5,'Events2',3);

bad = {};
for k = 1:numel(watch)
    nm = watch{k};
    w = which(nm);
    if isempty(w), continue, end
    if ~startsWith(lower(w),lower(folder))
        bad{end+1} = sprintf('    %-16s shadowed by %s',nm,w); %#ok<AGROW>
        continue
    end
    if isfield(arity,nm)
        try
            got = nargin(nm);
            if got ~= arity.(nm)
                bad{end+1} = sprintf(['    %-16s STALE LOADED COPY ' ...
                    '(takes %d args, project version takes %d)'], ...
                    nm,got,arity.(nm)); %#ok<AGROW>
            end
        catch
        end
    end
end
if ~isempty(bad)
    warning('setup_paths:shadowed', ...
        ['%d project file(s) are not the ones that will be used:\n%s\n' ...
         'Fix with:\n    clear classes; restoredefaultpath; setup_paths\n' ...
         '("clear" alone does NOT reload a classdef.)'], ...
        numel(bad), strjoin(bad,newline));
end
end
