classdef AeroMap
    %AEROMAP Interpolate CFD aero coefficients over F/R ride-height offsets.
    % CSV offsets are inches relative to the CFD reference ride heights.
    % evaluateNumeric provides a batched numeric interface for solver hot paths.

    properties (SetAccess = private)
        sourcePath
        frontRangeIn
        rearRangeIn
        interpolationMode = "scattered"
        claInterpolator
        cdaInterpolator
        copInterpolator
        claGridInterpolator
        cdaGridInterpolator
        copGridInterpolator
        sampleFrontOffsetIn
        sampleRearOffsetIn
        sampleCla
        sampleCda
        sampleCop
        claScale = 1
        cdaScale = 1
        copOffset = 0
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
            cla = cla(valid); cda = cda(valid); cop = cop(valid);
            obj.sourcePath = string(csvPath);
            obj.frontRangeIn = [min(front),max(front)];
            obj.rearRangeIn = [min(rear),max(rear)];
            obj.sampleFrontOffsetIn = front;
            obj.sampleRearOffsetIn = rear;
            obj.sampleCla = cla;
            obj.sampleCda = cda;
            obj.sampleCop = cop;

            % A complete rectilinear grid is safe for gridded interpolation
            % only when every cell is planar for all coefficients. In that
            % case bilinear interpolation matches the existing triangulated
            % linear scattered interpolation to floating-point precision.
            [useGrid,frontAxis,rearAxis,claGrid,cdaGrid,copGrid] = ...
                AeroMap.buildSafeRectangularGrid(front,rear,cla,cda,cop);
            if useGrid
                try
                    obj.claGridInterpolator = griddedInterpolant( ...
                        {frontAxis,rearAxis},claGrid,'linear','nearest');
                    obj.cdaGridInterpolator = griddedInterpolant( ...
                        {frontAxis,rearAxis},cdaGrid,'linear','nearest');
                    obj.copGridInterpolator = griddedInterpolant( ...
                        {frontAxis,rearAxis},copGrid,'linear','nearest');
                    obj.interpolationMode = "gridded";
                catch
                    % Keep the established scattered implementation as the
                    % safe fallback if MATLAB cannot construct a grid object.
                    obj.claGridInterpolator = [];
                    obj.cdaGridInterpolator = [];
                    obj.copGridInterpolator = [];
                end
            end

            % Nearest extrapolation makes a query outside the CFD envelope
            % explicit and bounded rather than inventing coefficients.
            obj.claInterpolator = scatteredInterpolant(front,rear,cla, ...
                'linear','nearest');
            obj.cdaInterpolator = scatteredInterpolant(front,rear,cda, ...
                'linear','nearest');
            obj.copInterpolator = scatteredInterpolant(front,rear,cop, ...
                'linear','nearest');
        end

        function aero = evaluate(obj,frontOffsetIn,rearOffsetIn)
            %EVALUATE Return the legacy coefficient struct for one state.
            [aero.cla,aero.cda,aero.D_f,aero.D_r,aero.outsideMap] = ...
                obj.evaluateNumeric(frontOffsetIn,rearOffsetIn);
            aero.frontOffsetIn = frontOffsetIn;
            aero.rearOffsetIn = rearOffsetIn;
        end

        function [cla,cda,D_f,D_r,outsideMap] = ...
                evaluateNumeric(obj,frontOffsetIn,rearOffsetIn)
            %EVALUATENUMERIC Evaluate one or more paired ride-height queries.
            % Inputs must have equal sizes, or one input may be scalar and
            % expand to the size of the other. Outputs preserve query shape:
            % CLA/CDA include configured scale factors, D_f/D_r are the
            % clamped aero load fractions, and outsideMap flags points beyond
            % the measured rectangular envelope. Out-of-envelope values use
            % bounded nearest-sample extrapolation, matching evaluate().
            [frontQuery,rearQuery] = AeroMap.expandQueryArrays( ...
                frontOffsetIn,rearOffsetIn);
            frontQuery = double(frontQuery);
            rearQuery = double(rearQuery);
            outsideMap = frontQuery < obj.frontRangeIn(1) | ...
                frontQuery > obj.frontRangeIn(2) | ...
                rearQuery < obj.rearRangeIn(1) | ...
                rearQuery > obj.rearRangeIn(2);

            claRaw = nan(size(frontQuery));
            cdaRaw = nan(size(frontQuery));
            copRaw = nan(size(frontQuery));
            finiteQuery = isfinite(frontQuery) & isfinite(rearQuery);

            if obj.interpolationMode == "gridded"
                inEnvelope = finiteQuery & ~outsideMap;
                if any(inEnvelope(:))
                    claRaw(inEnvelope) = obj.claGridInterpolator( ...
                        frontQuery(inEnvelope),rearQuery(inEnvelope));
                    cdaRaw(inEnvelope) = obj.cdaGridInterpolator( ...
                        frontQuery(inEnvelope),rearQuery(inEnvelope));
                    copRaw(inEnvelope) = obj.copGridInterpolator( ...
                        frontQuery(inEnvelope),rearQuery(inEnvelope));
                end

                outsideFinite = finiteQuery & outsideMap;
                for queryIndex = reshape(find(outsideFinite),1,[])
                    distanceSquared = ...
                        (obj.sampleFrontOffsetIn-frontQuery(queryIndex)).^2 + ...
                        (obj.sampleRearOffsetIn-rearQuery(queryIndex)).^2;
                    [~,nearestIndex] = min(distanceSquared);
                    claRaw(queryIndex) = obj.sampleCla(nearestIndex);
                    cdaRaw(queryIndex) = obj.sampleCda(nearestIndex);
                    copRaw(queryIndex) = obj.sampleCop(nearestIndex);
                end

                % Preserve the established scattered interpolant's behavior
                % for NaN/Inf inputs rather than inventing grid semantics.
                nonfiniteQuery = ~finiteQuery;
                if any(nonfiniteQuery(:))
                    claRaw(nonfiniteQuery) = obj.claInterpolator( ...
                        frontQuery(nonfiniteQuery),rearQuery(nonfiniteQuery));
                    cdaRaw(nonfiniteQuery) = obj.cdaInterpolator( ...
                        frontQuery(nonfiniteQuery),rearQuery(nonfiniteQuery));
                    copRaw(nonfiniteQuery) = obj.copInterpolator( ...
                        frontQuery(nonfiniteQuery),rearQuery(nonfiniteQuery));
                end
            else
                claRaw = obj.claInterpolator(frontQuery,rearQuery);
                cdaRaw = obj.cdaInterpolator(frontQuery,rearQuery);
                copRaw = obj.copInterpolator(frontQuery,rearQuery);
            end

            cla = obj.claScale*claRaw;
            cda = obj.cdaScale*cdaRaw;
            D_f = min(max(copRaw/100 + obj.copOffset,0),1);
            D_r = 1-D_f;
        end

        function obj = withCorrections(obj,claScale,cdaScale,copOffset)
            validateattributes(claScale,{'numeric'},{'scalar','finite','positive'}, ...
                mfilename,'claScale');
            validateattributes(cdaScale,{'numeric'},{'scalar','finite','positive'}, ...
                mfilename,'cdaScale');
            validateattributes(copOffset,{'numeric'},{'scalar','finite'}, ...
                mfilename,'copOffset');
            obj.claScale = claScale;
            obj.cdaScale = cdaScale;
            obj.copOffset = copOffset;
        end
    end

    methods (Static, Access = private)
        function [isSafe,frontAxis,rearAxis,claGrid,cdaGrid,copGrid] = ...
                buildSafeRectangularGrid(front,rear,cla,cda,cop)
            frontAxis = unique(front(:));
            rearAxis = unique(rear(:));
            frontCount = numel(frontAxis);
            rearCount = numel(rearAxis);
            isSafe = frontCount >= 2 && rearCount >= 2 && ...
                frontCount*rearCount == numel(front);
            claGrid = [];
            cdaGrid = [];
            copGrid = [];
            if ~isSafe, return, end

            [~,frontIndex] = ismember(front(:),frontAxis);
            [~,rearIndex] = ismember(rear(:),rearAxis);
            if any(frontIndex == 0 | rearIndex == 0) || ...
                    size(unique([frontIndex,rearIndex],'rows'),1) ~= numel(front)
                isSafe = false;
                return
            end

            gridSize = [frontCount,rearCount];
            linearIndex = sub2ind(gridSize,frontIndex,rearIndex);
            claGrid = nan(gridSize);
            cdaGrid = nan(gridSize);
            copGrid = nan(gridSize);
            claGrid(linearIndex) = cla(:);
            cdaGrid(linearIndex) = cda(:);
            copGrid(linearIndex) = cop(:);
            isSafe = AeroMap.isPlanarGrid(claGrid) && ...
                AeroMap.isPlanarGrid(cdaGrid) && AeroMap.isPlanarGrid(copGrid);
            if ~isSafe
                claGrid = [];
                cdaGrid = [];
                copGrid = [];
            end
        end

        function tf = isPlanarGrid(values)
            mixedDifference = values(1:end-1,1:end-1) + ...
                values(2:end,2:end) - values(2:end,1:end-1) - ...
                values(1:end-1,2:end);
            tolerance = 64*eps(max(1,max(abs(values(:)))));
            tf = all(abs(mixedDifference(:)) <= tolerance);
        end

        function [frontQuery,rearQuery] = expandQueryArrays(frontQuery,rearQuery)
            if ~isnumeric(frontQuery) || ~isnumeric(rearQuery)
                error('AeroMap:invalidQuery','Ride-height offsets must be numeric.');
            end
            if isscalar(frontQuery) && ~isscalar(rearQuery)
                frontQuery = repmat(frontQuery,size(rearQuery));
            elseif isscalar(rearQuery) && ~isscalar(frontQuery)
                rearQuery = repmat(rearQuery,size(frontQuery));
            elseif ~isequal(size(frontQuery),size(rearQuery))
                error('AeroMap:querySizeMismatch', ...
                    'Front and rear offsets must have equal sizes or be scalar.');
            end
        end
    end
end
