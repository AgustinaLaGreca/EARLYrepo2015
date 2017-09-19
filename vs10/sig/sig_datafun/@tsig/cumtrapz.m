function Y = cumtrapz(T);
% tsig/cumtrapz - CUMTRAPZ for tsig objects.
%    C=cumtrapz(S) for tsig object S is a tsig whose channel C(k) contains
%    an approximation of time integral of S's kth channel, S(k). C is
%    properly spaced, i.e. C =~ cumsum(S)*S.dt.
%  
%    To compute the definite time integral of T, see tsig/trapz.
% 
%    See also tsig/cumsum, tsig/dt, tsig/sum, tsig/plus.

if isvoid(T),
    error('Invalid input: void tsig object.'); 
end

% define portable cumtrapz function with proper scaling
DT = DT(T);
my_cumtrapz = @(x)cumtrapz(x*DT);
Y = wavefun(my_cumtrapz, T, 'tsig');



