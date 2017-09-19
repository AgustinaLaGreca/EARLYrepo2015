function tau = getWBdelay(S);
% transfer/getWBdelay - get wideband delay [ms] of transfer function
%    getWBdelay(S) returns the wideband delay in ms of Transfer object S.
%    If S is an array, an array of delays is returned.
%
%    See Transfer, Transfer/measure, Transfer/setWBdelay.

Sz = size(S);
for ii=1:numel(S),
    ttt = S(ii).WB_delay_ms;
    if isempty(ttt), ttt = nan; end
    tau(ii) = ttt;
end
tau = reshape(tau, Sz);
















