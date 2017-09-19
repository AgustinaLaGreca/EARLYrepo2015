function ir=isready(R,flag);
% playsound/isready - test if D/A is finished
%    isready(R) returns True of R is finished playing sound. If R has been
%    stopped, isready returns False: stopping is not the same as finishing.
%
%    isready(R,'sloppy') skips testing whether R's status is 'started' or 
%    'stopped'.
%
%    See also playsound/prepare, playsound/stop.

if nargin<2, flag=''; end

R = download(R);

if isequal('stopped', status(R)), % stopping ~= finishing; see help
    ir = 0;
    return;
end

if ~isequal('started', status(R)) && ~isequal('sloppy', flag),
    error('Can only query readiness of an object whose status is ''started''.');
end

St = seqplaystatus;
ir = ~St.Active;




