% =========================================================================
% -- Empirical CDF Generator
% -------------------------------------------------------------------------
% Description:
%   Evaluate the empirical cumulative distribution function of input
%   samples on a specified grid.
%
% Input:
%   x: Sampling grid for the empirical CDF.
%   in: Input samples used to form the empirical CDF.
%
% Output:
%   cdfdata: Empirical CDF values on the supplied grid.
%
% -- (c) 2008-2019 Christoph Studer
% =========================================================================

function cdfdata = cdfgen(x,in)
cdfdata = zeros(length(x),1);
    for i=1:length(x)
        cdfdata(i) = sum(in<x(i))/length(in);
    end
return
end

