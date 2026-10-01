function residual = qssConstraintResidual(c,ceq)
%QSSCONSTRAINTRESIDUAL Return a finite scalar residual for a QSS candidate.
inequalityValues = c(:);
equalityValues = ceq(:);
if isempty([inequalityValues;equalityValues]) || ...
        any(~isfinite([inequalityValues;equalityValues]))
    residual = Inf;
else
    % Inequalities are feasible on the negative side; only their positive
    % violation contributes. Equalities are signed residuals.
    residual = max([inequalityValues;abs(equalityValues);0]);
end
end
