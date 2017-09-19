function P = or(S,T);
% tsig/or - OR for tsig objects
%    S | T performs a sample-wise, channel-wise OR between tsig
%    objects S and T and returns a logical-valued tsig containing the
%    result of the AND operation. S and T must have the same number of
%    channels, unless one of them is single-channel. 
%
%    If the time domains of the corresponding channels of S and T are 
%    not identical, those parts that are missing in S but present in T are
%    provided, and set to zero (or false in the case if haslogic(S)). 
%    The same applies to those parts missing in T but present in S.
%    This convention is identical to the assumption of zeros (false)
%    outside the domain as used in tsig/plus.
%
%    S | X and X | S compare tsig S with logical X and return the result in a
%    logical tsig. X may be a scalar or a 1 x S.nchan array, in which case
%    the operation is performed channel-wise.
%
%    See also tsig/plus, tsig/relop.

op = eval(['@' mfilename]); 
% delegate the work to generic private relop
[P, Mess] = relop(S,T,op);
error(Mess);



