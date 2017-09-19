function [D, DS, S] = lastABF;
% lastABF - most recently loaded ABF data
%   [D, DS, s]  = lastABF()  
%   returns the ABF data D that was most recently loaded using readTKabf,
%   the corrsponding SGSRdataset DS and the TKpool entry s.
%
%    See also readTKABF, clampfit.

qq=MyFlag('LastRecordingLoaded');
[ExpID, RecID, icond] = deal(qq{:});

[D, DS, S] = readTKABF(ExpID, RecID, icond);

