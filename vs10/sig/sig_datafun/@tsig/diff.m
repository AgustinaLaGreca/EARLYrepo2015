function Y = diff(T);
% tsig/diff - DIFF for tsig objects: approximate time derivative.
%    D=diff(S) is the approximate time derivative of S, dS/dt, computed by
%    taking the difference of succesive samples and dividing them by S.dt.
%    By convention, D runs from S.onset to S.offset-S.dt (this is because 
%    onsets are always rounded to an integer number of sample periods).
% 
%    See also tsig/cumtrapz, tsig/onset.

if isvoid(T),
    error('Invalid input: void tsig object.'); 
end

DT = DT(T);
my_diff = @(x)diff(x*DT); % properly scaled version of diff
Y = wavefun(my_diff, T, 'tsig');



