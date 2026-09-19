function run = makeFixtureRun(caseInfo,varargin)
%MAKEFIXTURERUN Return a complete normalized run for injected-run tests.

if nargin < 1 || isempty(caseInfo)
    caseInfo = struct('id',"fixture",'label',"fixture",'carRole',"lap");
end
type = "lateral";
mode = "coast";
if nargin >= 2 && ~isempty(varargin{1})
    type = string(varargin{1});
end
if nargin >= 3 && ~isempty(varargin{2})
    mode = string(varargin{2});
end
fixture = makeRampFixture();
if type == "longitudinal"
    raw = fixture.legacyLongitudinalResult;
else
    raw = fixture.legacyRampResult;
end
settings = raw.settings;
run = rampSpeed.normalizeRampResult(raw,type, ...
    settings,caseInfo,struct('source',"makeFixtureRun", ...
    'fixture',true,'mode',mode));
run.mode = mode;
run.status = "completed";
fixedTimestamp = datetime(2026,1,1,0,0,0);
run.runMeta.created = fixedTimestamp;
run.runMeta.completed = fixedTimestamp;
run.runMeta.source = "makeFixtureRun";
end
