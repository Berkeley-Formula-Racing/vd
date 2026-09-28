function tests = test_rampSpeedAppAnalysisLayoutConfig
tests = functiontests(localfunctions);
end

function testAnalysisPagesUseOneCodeLevelLayoutDefinition(testCase)
root = fileparts(fileparts(mfilename("fullpath")));
archive = fullfile(root,"RampSpeedApp.mlapp");
temporaryRoot = tempname;
mkdir(temporaryRoot);
cleanup = onCleanup(@()rmdir(temporaryRoot,"s")); %#ok<NASGU>
unzip(archive,temporaryRoot);

documentPath = fullfile(temporaryRoot,"matlab","document.xml");
document = fileread(documentPath);

verifyNotEmpty(testCase,regexp(document, ...
    'function\s+layout\s*=\s*analysisPageLayout',"once"));
verifyNotEmpty(testCase,regexp(document, ...
    'function\s+axesArray\s*=\s*createAnalysisAxes',"once"));
verifyEmpty(testCase,regexp(document,'balanceMetricIds\s*=',"once"));
verifyEmpty(testCase,regexp(document,'suspensionMetricIds\s*=',"once"));
verifyEmpty(testCase,regexp(document,'balanceTitles\s*=',"once"));
verifyEmpty(testCase,regexp(document,'suspensionTitles\s*=',"once"));
end

function testBalanceAndSuspensionRemainTwoByTwo(testCase)
ensureRampSpeedAppPath();
app = RampSpeedApp("Visible","off");
cleanup = onCleanup(@()deleteIfValid(app)); %#ok<NASGU>

verifyEqual(testCase,numel(app.BalanceAxes),4);
verifyEqual(testCase,numel(app.SuspensionAxes),4);
verifyEqual(testCase,arrayfun(@(ax)ax.Layout.Row,app.BalanceAxes(:)), ...
    [1;1;2;2]);
verifyEqual(testCase,arrayfun(@(ax)ax.Layout.Column,app.BalanceAxes(:)), ...
    [1;2;1;2]);
verifyEqual(testCase,arrayfun(@(ax)ax.Layout.Row,app.SuspensionAxes(:)), ...
    [1;1;2;2]);
verifyEqual(testCase,arrayfun(@(ax)ax.Layout.Column,app.SuspensionAxes(:)), ...
    [1;2;1;2]);
end

function ensureRampSpeedAppPath()
root = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(root));
end

function deleteIfValid(app)
if isvalid(app)
    delete(app);
end
end
