function h2=has2varparams(D);
% Dataset/has2varparams - true when two stim params were varied independently
%   has2varparams(D) returns True if the stimulus presentation of dataset 
%    D prescribes the independent variation of two stimulus 
%    parameters, False otherwise. Arrays D allowed.
%
%   See also stimpresentx/has2varparams.

for ii=1:numel(D),
    h2(ii) = has2varparams(D(ii).Stim.Presentation);
end
h2 =  reshape(h2, size(D));


