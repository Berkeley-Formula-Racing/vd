function opts = setOptimoptions(numIter,constraintTol)
% inputs  numIter:       MaxFunctionEvaluations
%         constraintTol: ConstraintTolerance (optional, default 1e-2)
%
% 1e-2 is loose for constraints that include wheel torque residuals in N-m:
% a wide set of states qualifies as feasible, so two solves from slightly
% different starts can settle on visibly different answers. That is tolerable
% for a single g-g point, which is read on its own, and NOT tolerable inside
% vel_cornering_sweep, where 556 solves are chained and the spread compounds.
% Callers that need a repeatable answer pass a tighter value.

if nargin < 2 || isempty(constraintTol), constraintTol = 1e-2; end

% default algorithm is interior-point
opts = optimoptions('fmincon','MaxFunctionEvaluations',numIter, ...
    'ConstraintTolerance',constraintTol,'StepTolerance',1e-10,'Display','off');


% exitflag meaning: 1 = converged, 2 = change in x less than step tolerance
%   (optimality condition not fulfilled, but solution still found
%   0 = function evaluations exceeded (not converging)
%   -2 = no feasible point found 