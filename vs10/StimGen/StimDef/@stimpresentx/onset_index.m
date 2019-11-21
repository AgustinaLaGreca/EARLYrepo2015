function iSam = onset_index(SPx)
% stimPresentx/onset_index - sample indices at start of stim presentations
%    onset_index(SP) returns an array holding the zero based sample indices
%    of all the stimulus presentations, baselines included (first & last elements). 
%    This is the sample count ("tick value") at the start of a
%    presentation, which may or may not coincide with the onset of the
%    sound presentation, depending on the 'delay' value specified in the
%    stimulus menu.

iSam = SPx.SamOffset(1:end-1);







