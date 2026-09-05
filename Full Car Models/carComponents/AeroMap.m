classdef AeroMap
    %AEROMAP Interpolate CFD aero coefficients over F/R ride-height offsets.
    % CSV offsets are inches relative to the CFD reference ride heights.

    properties (SetAccess = private)
        sourcePath
        frontRangeIn
        rearRangeIn
        claInterpolator
        cdaInterpolator
        copInterpolator
    end

    methods
        function obj = AeroMap(csvPath)
            if ~isfile(csvPath)
                error('AeroMap:fileNotFound','Aeromap file not found: %s',csvPath);
            end
            T = readtable(csvPath,'VariableNamingRule','preserve');
            required = ["FFR offset","RRH offset","CLA","CDA","COP"];
            if ~all(ismember(required,string(T.Properties.VariableNames)))
                error('AeroMap:missingColumns', ...
                    'Aeromap must contain FFR offset, RRH offset, CLA, CDA, and COP columns.');
            end

            front = double(T.("FFR offset"));
            rear  = double(T.("RRH offset"));
            cla   = double(T.CLA);
            cda   = double(T.CDA);
            cop   = double(T.COP);
            valid = isfinite(front) & isfinite(rear) & isfinite(cla) & ...
                isfinite(cda) & isfinite(cop);
            if nnz(valid) < 3
                error('AeroMap:notEnoughPoints','Aeromap needs at least three finite samples.');
            end

            front = front(valid); rear = rear(valid);
            obj.sourcePath = string(csvPath);
            obj.frontRangeIn = [min(front),max(front)];
            obj.rearRangeIn = [min(rear),max(rear)];
            % Nearest extrapolation makes a query outside the CFD envelope
            % explicit and bounded rather than inventing coefficients.
            obj.claInterpolator = scatteredInterpolant(front,rear,cla(valid), ...
                'linear','nearest');
            obj.cdaInterpolator = scatteredInterpolant(front,rear,cda(valid), ...
                'linear','nearest');
            obj.copInterpolator = scatteredInterpolant(front,rear,cop(valid), ...
                'linear','nearest');
        end

        function aero = evaluate(obj,frontOffsetIn,rearOffsetIn)
            aero.cla = obj.claInterpolator(frontOffsetIn,rearOffsetIn);
            aero.cda = obj.cdaInterpolator(frontOffsetIn,rearOffsetIn);
            aero.D_f = min(max(obj.copInterpolator(frontOffsetIn,rearOffsetIn)/100,0),1);
            aero.D_r = 1-aero.D_f;
            aero.frontOffsetIn = frontOffsetIn;
            aero.rearOffsetIn = rearOffsetIn;
            aero.outsideMap = frontOffsetIn < obj.frontRangeIn(1) || ...
                frontOffsetIn > obj.frontRangeIn(2) || ...
                rearOffsetIn < obj.rearRangeIn(1) || rearOffsetIn > obj.rearRangeIn(2);
        end
    end
end
