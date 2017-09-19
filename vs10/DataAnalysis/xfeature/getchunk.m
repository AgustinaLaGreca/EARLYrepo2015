function s = getchunk(S,C);
% Getchunk - get single segment from signal
%   Getchunk(S,C) returns the segment S(C.istart:C.iend).
%   C must be single struct with fields istart & iend.
%   If S is a matrix, the columns are interpreted as time series and
%   selection is done in the 1st dimension.
%
%   Helper function for EvalChuncks.
%
%   Seel also EvalChuncks.

if isvector(S),
    s = S(C.istart:C.iend);
else,
    s = S(C.istart:C.iend,:);
end

