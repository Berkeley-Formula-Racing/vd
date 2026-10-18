function [setups,designTable] = sampleSetups(baseSetup,samplingType,samples)
%SAMPLESETUPS Sample editable setups using the legacy QSS table semantics.
if nargin < 2 || isempty(samplingType), samplingType = "FullFactorial"; end
if nargin < 3, samples = []; end
if ~isstruct(baseSetup) || ~isscalar(baseSetup) || ...
        ~isfield(baseSetup,'parameters') || ~isstruct(baseSetup.parameters)
    error('vehicleConfig:badBaseSetup', ...
        'baseSetup must be a scalar setup struct with a parameters struct.');
end

names = vehicleConfig.parameterNames();
values = cell(1,numel(names));
for j = 1:numel(names)
    values{j} = baseSetup.parameters.(names{j});
end
nVar = numel(names);
switch lower(string(samplingType))
    case "explicit"
        sampleInput = samples;
        if ~istable(sampleInput)
            error('parameters_loop:badExplicitInput', ...
                'Explicit sampling requires a table of parameter overrides.');
        end
        if any(~cellfun(@isscalar,values))
            error('parameters_loop:nonScalarBaseline', ...
                'Explicit sampling requires scalar baseline values.');
        end
        P = repmat(cellfun(@(x) x(1),values),height(sampleInput),1);
        for k = 1:width(sampleInput)
            j = find(strcmp(names,sampleInput.Properties.VariableNames{k}),1);
            if isempty(j)
                error('parameters_loop:unknownExplicitParameter', ...
                    'Unknown explicit parameter %s.', ...
                    sampleInput.Properties.VariableNames{k});
            end
            P(:,j) = sampleInput{:,k};
        end

    case {"fullfactorial","grid"}
        grid = cell(1,nVar);
        [grid{:}] = ndgrid(values{:});
        P = zeros(numel(grid{1}),nVar);
        for j = 1:nVar
            P(:,j) = grid{j}(:);
        end

    case {"lhs","random"}
        validateattributes(samples,{'numeric'},{'scalar','integer','positive'});
        U = rand(samples,nVar);
        if strcmpi(samplingType,"LHS")
            for j = 1:nVar
                U(:,j) = (randperm(samples)' - U(:,j))/samples;
            end
        end
        P = zeros(samples,nVar);
        for j = 1:nVar
            bounds = values{j};
            if isscalar(bounds)
                P(:,j) = bounds;
            else
                lo = min(bounds);
                hi = max(bounds);
                P(:,j) = lo + U(:,j)*(hi-lo);
            end
        end

    otherwise
        error('parameters_loop:badSamplingType', ...
            'samplingType must be FullFactorial, LHS, or Random.');
end

designTable = array2table(P,'VariableNames',names);
setups = repmat(baseSetup,size(P,1),1);
for i = 1:size(P,1)
    for j = 1:nVar
        setups(i).parameters.(names{j}) = P(i,j);
    end
    setupId = sprintf('case-%06d',i);
    setups(i).id = string(setupId);
    setups(i).label = string(sprintf('Case %06d',i));
end
end
