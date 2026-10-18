function [lltd] = LLTDVals(m,muf,mur,msf,msr,hg,huf,hur,zf,zr,t,l,kc,trs,rsd)
%LLTDVALS Summary of this function goes here
%   Detailed explanation goes here
lambda = rsd;
mu = kc/trs;

lltd = ((lambda^2 - (mu + 1) * lambda) * dsf * msf) / ((lambda^2 - lambda - mu) * hg * m) - ...
             (mu * lambda * dsr * msr) / ((lambda^2 - lambda - mu) * hg * m) + ...
             zf * msf / (hg * m) + huf * muf / (hg * m);
end

