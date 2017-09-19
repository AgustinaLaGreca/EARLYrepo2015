function N = Ncond(DS);
% Dataset/Ncond - # stimulus conditions in DS, whether recorded or not
%    Ncond(DS) returns the number of stimulus conditions specified at
%    recording time. Not that interrupted recordings may not hold data for
%    all stimulus conditions.
%
%    See also StimPresent, Dataset/NrepRec.

N = DS.Stim.Presentation.Ncond;







