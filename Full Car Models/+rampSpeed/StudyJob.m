classdef StudyJob < handle
    %STUDYJOB Mutable lifecycle state shared by the study executor and UI.

    properties
        state = "idle"
        study = struct()
        progressEvents = struct.empty(0,1)
        error = []
        cancelRequested = false
        future = []
        cancelFile = ""
        started = NaT
        completed = NaT
    end

    methods
        function obj = StudyJob(cancelFile)
            if nargin >= 1 && ~isempty(cancelFile)
                obj.cancelFile = string(cancelFile);
            end
            obj.started = datetime("now");
        end

        function requestCancel(obj)
            obj.cancelRequested = true;
            if strlength(string(obj.cancelFile)) == 0
                return
            end
            try
                fid = fopen(char(obj.cancelFile),"w");
                if fid >= 0
                    fclose(fid);
                end
            catch
                % The in-memory token remains authoritative for serial runs.
            end
        end

        function value = isCancelled(obj)
            value = logical(obj.cancelRequested);
            if value || strlength(string(obj.cancelFile)) == 0
                return
            end
            try
                value = isfile(char(obj.cancelFile));
            catch
                value = false;
            end
        end

        function delete(obj)
            fileName = string(obj.cancelFile);
            if strlength(fileName) > 0 && isfile(char(fileName))
                try
                    delete(char(fileName));
                catch
                end
            end
        end
    end
end
