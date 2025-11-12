function [ F,J] = gausslinefun( params, x)
%DOGFUN evaluates modified gaussian with params at x
%   Output: 
%       F: function value
%       J: Jacobian of the parameters
%Written by Dimos.
%==========================================================================
%Unwrap inputs
%==========================================================================
x = x(:);

mx=params(1); ss = params(2); Ac=params(3); 
%==========================================================================
%Calculate function value
%==========================================================================
inExpC = -(x - mx).^2/(2 * ss^2);
F      = Ac * exp(inExpC);
%==========================================================================
%Calculate Jacobian
%==========================================================================
if nargout>1
    
    dmx = (x-mx)/ss^2;
    J1  = Ac*exp(inExpC).*dmx;

    dss = ss^-3*(x-mx).^2;
    J2  = Ac*exp(inExpC).*dss ;
    J3 = exp(inExpC); 
    
    J=[J1 J2 J3];
end
%==========================================================================
end

