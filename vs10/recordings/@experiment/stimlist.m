function [L, R, T] = stimlist(E, Ename);
% experiment/stimlist - list stimulus conditions of all recordings of experiment
%    stimlist(E) or stimlist(experiment, 'Expname') returns a struct array
%    whose elements are the stimlists of all the datasets in experiment E.
%
%    [L, R, T] = stimlist(...) also returns elementary recording info in struct
%    array R and the stimulus types in cell string T.
%
%    See also dataset/stimlist.

if isvoid(E),
    E = dummy(E, Ename);
end
Stat = status(E);
iRec = [Stat.AllSaved.iDataset]; % indices of all datasets
DS = dummy(dataset, name(E), iRec);
[L,R,T] = stimlist(DS);












