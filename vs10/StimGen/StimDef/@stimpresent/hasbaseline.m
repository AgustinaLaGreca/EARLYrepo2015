function h2 = hasbaseline(P)
% stimPresent/hasbaseline - true when stimulus includes baselines
%    hasbaseline(P) returns True if stimpresent object P includes baseline
%    measurements.

for ii=1:numel(P),
    h2(ii) = isfield(struct(P(ii)), 'PreDur');
end
h2 = reshape(h2,size(P));


