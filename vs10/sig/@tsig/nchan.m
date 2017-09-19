function N = nchan(T)
% tsig/nchan - number of channels of tsig object.
%    nchan(T) of T.nchan returns the number of chabbels in tsig object T.
%    By convention, nchan returns NaN for void tsig objects.
%
%    See also tsig, tsig/nchan, tsig/nchan, tsig/isvoid, "methods tsig".

if isvoid(T), N = nan;
else,
    N = size(T.t0, 2);
end




