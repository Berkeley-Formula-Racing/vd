function copy = duplicateSetup(source,newId,newLabel)
%DUPLICATESETUP Copy a serializable setup specification for editing.

if ~isstruct(source) || ~isscalar(source)
    error('rampSpeed:invalidSetup', ...
        'The source setup must be a scalar struct.');
end
copy = source;
if nargin < 2 || isempty(newId)
    newId = string(source.id) + "-copy";
end
if nargin < 3 || isempty(newLabel)
    newLabel = string(source.label) + " copy";
end
copy.id = string(newId);
copy.label = string(newLabel);
copy.source = "duplicate";
copy.isBaseline = false;
end
