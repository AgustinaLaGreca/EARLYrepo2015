function S = repmat(S,N,M);
% tsig/repmat - tile tsig objects.
%   R=repmat(S,M,N) uses tsig/horzcat and tsig/vertcat to concatenate the
%   channels of S N-fold, and to append the channels of S M-fold,
%   repectively. Note that The result R is still a *single* tsig object;
%   tsig arrays are not allowed.
%
%   If you do need arrays of tsog objects, use cell arrays.
%
%   See also tsig/subsref, tsig/vertcat, tsig/horzcat.

error(numericTest(N,'scalar/posint/real/noninf','Argument N is '))
error(numericTest(M,'scalar/posint/real/noninf','Argument M is '))

if N>1,
    SN = repmat({S},1,N);
    S = vertcat(SN{:});
end

if M>1,
    SM = repmat({S}, 1,M);
    S = horzcat(SM{:});
end


