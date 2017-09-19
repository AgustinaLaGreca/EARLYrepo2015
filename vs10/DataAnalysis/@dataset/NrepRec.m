function N=NrepRec(D);
% Dataset/NrepRec - # repetitions for each recorded stimulus condition
%   N = NrepRec(D) returns an array N whose k-th element N(k) holds the 
%   number of repetitions of for stimulus condition k for which recordings
%   are present in dataset D. If D is complete, N consists of all Nrep
%   values, where Nrep is the number of specified reps. If D is incomplete,
%   the number of res varies in a way that depends on the time of
%   interruption and the order of presentation.
%
%   See also dataset/NpresRecorded, dataset/iscomplete.


AllReps = D.Stim.Presentation.iCond(1:NpresRecorded(D));
N = hist(AllReps, 1:D.Stim.Presentation.Ncond+2); % allow for baselines
N = N(1:D.Stim.Presentation.Ncond);

if isfield(D.Stim, 'Stutter')
   if strcmpi(D.Stim.Stutter,'on')
       N(1) = N(1)-1;
   end
end


