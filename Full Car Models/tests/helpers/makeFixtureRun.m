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
settings = fixture.legacyRampResult.settings;
run = rampSpeed.normalizeRampResult(fixture.legacyRampResult,type, ...
    settings,caseInfo,struct('source',"makeFixtureRun", ...
    'fixture',true,'mode',mode));
run.mode = mode;
run.status = "completed";
run.runMeta.completed = datetime('now');
run.runMeta.source = "makeFixtureRun";
end
