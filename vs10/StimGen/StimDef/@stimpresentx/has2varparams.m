function h2 = has2varparams(P)
% stimPresentx/has2varparams - true when two stim params are varied 
%    has2varparams(P) returns True if stimpresentx object P prescribes the
%    independent variation of two stimulus parameters, False otherwise. 
%    This amounts to checking whether P.Y is non-empty. Arrays P allowed.
%
%    See also Dataset/has2varparams.

for ii=1:numel(P),
    h2(ii) = ~isempty(P(ii).Y);
end
h2 = reshape(h2,size(P));


