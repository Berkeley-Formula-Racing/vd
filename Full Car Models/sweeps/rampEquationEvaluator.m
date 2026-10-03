function evaluator = rampEquationEvaluator(car,useRideHeightContext)
% Cache only the equation state required by ramp constraints.
%
% Coast's reduced solve carries the coupled ride-height context between
% compatible states. Balanced retains the historical nine-output equation
% call in its full-state hot loop: constructing and validating a context for
% every full-state finite-difference probe costs more than the exact solve.
% Both paths memoize the exact state and keep full post-solve residual checks.

if nargin < 2 || isempty(useRideHeightContext)
    useRideHeightContext = true;
end
useRideHeightContext = logical(useRideHeightContext);

lastP = [];
lastState = struct();
lastRideHeightContext = [];
equationCallCount = 0;

evaluator = struct('evaluate',@evaluate,'equationCalls',@equationCalls);

    function state = evaluate(P)
        P = P(:).';
        if isempty(lastP) || ~isequal(P,lastP)
            if useRideHeightContext
                [engineRpm,beta,latAccel,longAccel,yawAccel,wheelAccel,~,~, ...
                    Fzvirtual,~,~,~,~,~,~,~,lastRideHeightContext] = ...
                    car.equations(P,lastRideHeightContext, ...
                        struct('qssValidityOnly',true));
            else
                % Request the same nine outputs as the historical balanced
                % ramp solver. This avoids constructing a ride-height context
                % and diagnostic tuple for every full-state probe.
                [engineRpm,beta,latAccel,longAccel,yawAccel,wheelAccel,~,~, ...
                    Fzvirtual] = car.equations(P);
                lastRideHeightContext = [];
            end
            lastP = P;
            equationCallCount = equationCallCount + 1;
            lastState = struct('engineRpm',engineRpm,'beta',beta, ...
                'latAccel',latAccel,'longAccel',longAccel, ...
                'yawAccel',yawAccel,'wheelAccel',wheelAccel, ...
                'Fzvirtual',Fzvirtual, ...
                'rideHeightContext',lastRideHeightContext);
        end
        state = lastState;
    end

    function count = equationCalls()
        count = equationCallCount;
    end
end
