function evaluator = steadyStateEvaluator(car)
%STEADYSTATEEVALUATOR Cache one exact vehicle-equation evaluation.
%   FMINCON commonly asks the objective and nonlinear constraints about the
%   same state in sequence. Reusing that exact result avoids solving the
%   coupled ride-height/aeromap system twice without approximating it.

lastP = [];
lastState = struct();
equationCallCount = 0;

evaluator = struct('evaluate',@evaluate,'equationCalls',@equationCalls);

    function state = evaluate(P)
        P = P(:).';
        if isempty(lastP) || ~isequal(P,lastP)
            [engineRpm,beta,latAccel,longAccel,yawAccel,wheelAccel,~,~,Fzvirtual] = ...
                car.equations(P);
            lastP = P;
            equationCallCount = equationCallCount + 1;
            lastState = struct('engineRpm',engineRpm,'beta',beta, ...
                'latAccel',latAccel,'longAccel',longAccel, ...
                'yawAccel',yawAccel,'wheelAccel',wheelAccel, ...
                'Fzvirtual',Fzvirtual);
        end
        state = lastState;
    end

    function count = equationCalls()
        count = equationCallCount;
    end
end
