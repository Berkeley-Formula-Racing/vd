classdef simLog
% Append-only run log for the lap sim. One row per top-level run -- a sweep, a
% full events pass, a standalone g-g -- recording how long it took and the
% context needed to read that number later: what was run, how many cases, how
% many workers, the resulting event times, the car, and the git state of the
% code that produced it.
%
% USAGE
%   job = simLog.start('grip sweep');
%   ... run ...
%   simLog.finish(job, 'events',{'skidpad','autocross'}, 'workers',8, ...
%                      'nCases',9, 'times',ev.times, 'car',car);
%
% start() stamps the clock and starts a timer; finish() stops it and appends
% one CSV row. Everything after the job handle is optional name/value context.
% Unrecognised names are not dropped -- they land in a free-text `details`
% column -- so logging too much never loses information.
%
%   simLog.show           print the last 20 rows
%   simLog.show(n)        print the last n
%   T = simLog.read       the whole log as a table
%   simLog.path           where the file lives (sim_log.csv, next to carConfig)
%
% ONE WRITER. Call this from the CLIENT, once, after a run finishes -- never
% from inside a parfor, where concurrent workers would interleave and corrupt
% the file. Two MATLAB sessions writing at once is also best avoided; each row
% is a single fprintf, so the worst a collision does is tear one line.

    properties (Constant, Access = private)
        COLS = {'timestamp','label','events','n_cases','workers','elapsed_s', ...
                's_per_case','skidpad','accel','autocross','endurance', ...
                'gripF','gripR','mass','git','dirty','host','matlab','details'};
    end

    methods (Static)
        function job = start(label)
            % Begin timing a run. LABEL names the kind of run for the log.
            if nargin < 1, label = 'run'; end
            job = struct('t0',tic,'started',datetime('now'),'label',char(label));
        end

        function row = finish(job, varargin)
            % Stop the timer started by start() and append a row. Optional
            % name/value context, all case-insensitive:
            %   events   cellstr or char   -> joined with '+'
            %   nCases   numeric           -> also drives s_per_case
            %   workers  numeric
            %   times    struct            -> pulls skidpad/accel/autocross/endurance
            %   car      Car               -> pulls mass, gripScaleF, gripScaleR
            %   skidpad/accel/autocross/endurance/gripF/gripR/mass  numeric overrides
            %   notes | details            text
            % anything else is appended to `details` as name=value.
            elapsed = toc(job.t0);

            d = struct('events','','n_cases',[],'workers',[],'skidpad',[], ...
                'accel',[],'autocross',[],'endurance',[],'gripF',[],'gripR',[], ...
                'mass',[],'details','');
            extra = {};
            k = 1;
            while k < numel(varargin)
                key = lower(char(varargin{k}));
                val = varargin{k+1};
                k = k + 2;
                switch key
                    case 'events'
                        if iscell(val), d.events = strjoin(cellfun(@char,val,'uni',0),'+');
                        else,           d.events = char(string(val)); end
                    case {'ncases','n_cases'}, d.n_cases = val;
                    case 'workers',            d.workers = val;
                    case 'times'
                        for fn = {'skidpad','accel','autocross','endurance'}
                            if isstruct(val) && isfield(val,fn{1}) && ~isempty(val.(fn{1}))
                                d.(fn{1}) = val.(fn{1})(1);
                            end
                        end
                    case 'car'
                        try
                            d.mass  = val.M;
                            d.gripF = val.gripScaleF;
                            d.gripR = val.gripScaleR;
                        catch
                        end
                    case 'skidpad',   d.skidpad   = val;
                    case 'accel',     d.accel     = val;
                    case 'autocross', d.autocross = val;
                    case 'endurance', d.endurance = val;
                    case 'gripf',     d.gripF     = val;
                    case 'gripr',     d.gripR     = val;
                    case 'mass',      d.mass      = val;
                    case {'notes','details'}, d.details = char(string(val));
                    otherwise
                        extra{end+1} = sprintf('%s=%s',key,simLog.scalarStr(val)); %#ok<AGROW>
                end
            end
            if ~isempty(extra)
                if isempty(d.details), d.details = strjoin(extra,'; ');
                else, d.details = [d.details '; ' strjoin(extra,'; ')]; end
            end

            if ~isempty(d.n_cases) && d.n_cases > 0, spc = elapsed/d.n_cases;
            else,                                    spc = elapsed; end

            [commit,dirty] = simLog.gitInfo();
            host = getenv('COMPUTERNAME');
            if isempty(host), host = getenv('HOSTNAME'); end
            if isempty(host), host = 'unknown'; end

            row = { ...
                char(string(job.started,'yyyy-MM-dd HH:mm:ss')), ...
                job.label, d.events, ...
                simLog.numStr(d.n_cases,'%d'), simLog.numStr(d.workers,'%d'), ...
                sprintf('%.2f',elapsed), sprintf('%.2f',spc), ...
                simLog.numStr(d.skidpad,'%.4f'), simLog.numStr(d.accel,'%.4f'), ...
                simLog.numStr(d.autocross,'%.4f'), simLog.numStr(d.endurance,'%.3f'), ...
                simLog.numStr(d.gripF,'%.4g'), simLog.numStr(d.gripR,'%.4g'), ...
                simLog.numStr(d.mass,'%.1f'), commit, mat2str(dirty), ...
                host, version('-release'), d.details};

            simLog.appendRow(row);
            if nargout == 0
                fprintf('  logged: %s  %.1f s  (%s) -> %s\n', ...
                    job.label, elapsed, d.events, simLog.path);
                clear row
            end
        end

        function p = path()
            % sim_log.csv sits next to carConfig (one up from utilities/).
            here = fileparts(mfilename('fullpath'));
            p = fullfile(fileparts(here),'sim_log.csv');
        end

        function T = read()
            % The whole log as a table, or an empty table if nothing logged.
            p = simLog.path;
            if ~isfile(p), T = table(); return, end
            opts = detectImportOptions(p,'Delimiter',',');
            % Force the text columns to string. Otherwise a short git hash like
            % 97e1622 is read as scientific notation (97e1622 = Inf), and a
            % details field can be mis-typed too. Numeric columns stay numeric.
            textCols = {'timestamp','label','events','git','dirty','host','matlab','details'};
            tc = intersect(textCols,opts.VariableNames);
            if ~isempty(tc), opts = setvartype(opts,tc,'string'); end
            T = readtable(p,opts);
        end

        function show(n)
            % Print the last n rows (default 20).
            if nargin < 1, n = 20; end
            T = simLog.read();
            if isempty(T), fprintf('sim log is empty (%s)\n',simLog.path); return, end
            fprintf('sim log: %d runs, %s\n',height(T),simLog.path);
            disp(T(max(1,end-n+1):end,:));
        end
    end

    methods (Static, Access = private)
        function appendRow(row)
            p = simLog.path;
            newfile = ~isfile(p);
            fid = fopen(p,'a');
            if fid < 0
                warning('simLog:cannotWrite','could not open %s for the run log',p);
                return
            end
            closer = onCleanup(@() fclose(fid));
            if newfile
                fprintf(fid,'%s\n',strjoin(simLog.COLS,','));
            end
            esc = cellfun(@simLog.csvEsc,row,'uni',0);
            fprintf(fid,'%s\n',strjoin(esc,','));
        end

        function s = csvEsc(v)
            s = char(string(v));
            if contains(s,{',','"',newline,char(13)})
                s = ['"' strrep(s,'"','""') '"'];
            end
        end

        function s = numStr(v,fmt)
            if isempty(v) || (isnumeric(v) && ~isfinite(v)), s = ''; else, s = sprintf(fmt,v); end
        end

        function s = scalarStr(v)
            if isnumeric(v) && isscalar(v), s = num2str(v);
            elseif ischar(v) || isstring(v), s = char(string(v));
            elseif islogical(v) && isscalar(v), s = mat2str(v);
            else, s = ['<' class(v) '>']; end
        end

        function [commit,dirty] = gitInfo()
            % Best effort. The commit is cached for the session (HEAD rarely
            % moves mid-run); dirty is re-checked each time, since editing code
            % between runs is exactly what makes a log row's provenance matter.
            persistent cachedCommit repoDir
            commit = ''; dirty = false;
            if ispc, dev = ' 2>NUL'; else, dev = ' 2>/dev/null'; end
            try
                if isempty(repoDir)
                    repoDir = fileparts(fileparts(mfilename('fullpath')));
                end
                if isempty(cachedCommit)
                    [st,out] = system(sprintf('git -C "%s" rev-parse --short HEAD%s',repoDir,dev));
                    if st == 0, cachedCommit = strtrim(out); else, cachedCommit = ''; end
                end
                commit = cachedCommit;
                [st2,out2] = system(sprintf('git -C "%s" status --porcelain%s',repoDir,dev));
                dirty = (st2 == 0) && ~isempty(strtrim(out2));
            catch
            end
        end
    end
end
