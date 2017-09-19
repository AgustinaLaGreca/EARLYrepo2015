function T = JLtimeSinceRecstart(Urecidx);
% JLtimeSinceRecstart - time in minutes since start of recording
%   JLtimeSinceRecstart(Urecidx)
%
%   See also JLrecStartTime

% try cached value
[T, CFN, CP] = getcache(mfilename, Urecidx);
if ~isempty(T), return; end

if numel(Urecidx)>1, % recursive handling
    for ii=1:numel(numel(Urecidx)),
        T(ii) = JLtimeSinceRecstart(Urecidx(ii));
    end
    putcache(CFN, 5e3, CP, T);
    return;
end


% find it from data
D = JLdatastruct(Urecidx);


T = D.abfdatenum - JLrecStartTime(D.icell_run); % age in days
T = T*24*60; % age in minutes

putcache(CFN, 5e3, CP, T);





