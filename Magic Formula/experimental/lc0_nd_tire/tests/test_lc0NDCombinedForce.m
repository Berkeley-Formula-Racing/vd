function tests = test_lc0NDCombinedForce
tests = functiontests(localfunctions);
end

function testPreservesPureAxisAndCapsCombinedDemand(testCase)
% Break caught: a combined-slip cap must not degrade a feasible pure force,
% and an overloaded combined demand must remain inside its stated ellipse.
addpath(fileparts(fileparts(mfilename('fullpath'))));
Fz = 100;
muX = 1.2;
muY = 1.0;
p = 2;

[fxPure,fyPure,uPure] = lc0NDCombinedForce(60,0,muX,muY,Fz,p);
verifyEqual(testCase,fxPure,60,'AbsTol',1e-12);
verifyEqual(testCase,fyPure,0,'AbsTol',1e-12);
verifyEqual(testCase,uPure,0.5,'AbsTol',1e-12);

[fx,fy,u] = lc0NDCombinedForce(120,100,muX,muY,Fz,p);
verifyGreaterThan(testCase,u,1);
verifyLessThanOrEqual(testCase, ...
    hypot(fx/(muX*Fz),fy/(muY*Fz)),1+1e-12);
verifyGreaterThan(testCase,fx,0);
verifyGreaterThan(testCase,fy,0);
end
