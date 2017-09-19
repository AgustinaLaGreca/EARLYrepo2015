function h2 = hasbaseline(P)
% stimPresentx/hasbaseline - true when stimulus includes baselines
%    hasbaseline(P) returns True if stimpresentx object includes baseline
%    measurements.

for ii=1:numel(P),
    h2(ii) = isfield(struct(P(ii)), 'PreDur');
end
h2 = reshape(h2,size(P));


