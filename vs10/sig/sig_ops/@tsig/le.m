function P = le(S,T);
% tsig/le - LE for tsig objects
%    S<=T performs a sample-wise, channel-wise comparison between tsig
%    objects S and T and returns a logical-valued tsig containing the
%    result of the comparison. S and T must have the same number of
%    channels, unless one of them is single-channel. 
%
%    If the time domains of the corresponding channels of S and T are 
%    not identical, those parts that are missing in S but present in T are
%    provided, and set to zero (or false in the case if haslogic(S)). 
%    The same applies to those parts missing in T but present in S.
%    This convention is identical to the assumption of zeros outside the 
%    domain as used in tsig/plus.
%
%    S<=X and X<=S compare tsig S with numerical X and return the result in a
%    logical tsig. X may be a scalar or a 1 x S.nchan array, in which case
%    the comparison is done channel-wise.
%
%    See also tsig/plus.

op = eval(['@' mfilename]); 
% delegate the work to generic private relop
[P, Mess] = relop(S,T,op);
error(Mess);



