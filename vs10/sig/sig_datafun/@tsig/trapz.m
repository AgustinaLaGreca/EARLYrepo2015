function Y = trapz(T,t0,t1);
% tsig/trapz - TRAPZ for tsig objects.
%    trapz(S) computes an approximation of the time integral of tsig S
%    via the trapezoidal method with proper time spacing in ms. That is, 
%    C =~ sum(S)*S.dt. The integral is computed independently in each 
%    channel of S and returned in a row vector C. In each channel the time
%    range of the integral is the channel's own time range.
% 
%    Trapz(S,t0,t1) restricts the integral to the time interval t0 to t1 ms. 
%    Zeros are padded for any time intervals outside a channels time range.
%    t0 and t1 may be 1 x nchan(S) row vectors. See tsig/cut for special 
%    values of t0 and t1.
% 
%    See also tsig/cumsum, tsig/cut, tsig/dt.

if isvoid(T),
    error('Invalid input: void tsig object.'); 
end

if nargin>1, % restrict the range
    try, 
        T = cut(T,t0,t1); % this will automatically pad zeros where needed
    catch, 
        error(lasterr);
    end
end

DT = DT(T);
my_trapz = @(x)trapz(x*DT); % properly scaled version of trap
Y = wavefun(my_trapz, T, 'cat');



