function Nrep = RueNrep(FN);
% RueNrep - rep count of Rue recording
%    Nrep = RueNrep('p90305Bc') returns the # reps of recording p90305Bc.
%
%    See also RueList.

if ~isequal('c', lower(FN(end))) && ~isequal('i', lower(FN(end))),
    FN = [FN 'i'];
end
   
[RL, Nrep] = RueList;
if ~ismember(FN, RL),
    error(['Cannot find data for cell ''' FN '''.']);
end
imatch = strmatch(FN, RL);
Nrep = Nrep(imatch);






