function [c,ceq] = steadyStateConstraint4(state,P,latAccelTarget)
%STEADYSTATECONSTRAINT4 Constraint4 from a cached vehicle-equation state.

c = [state.engineRpm-13000,abs(state.beta)-20,-state.Fzvirtual(1:4)];
ceq = [P(3)*P(5)-latAccelTarget,state.latAccel,state.yawAccel, ...
    state.wheelAccel(1:4)];
end
