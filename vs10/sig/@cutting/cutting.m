function U=cutting(X)
% cuttin- cutting constructor
%    For internal use only. Called by tsig/subsref
%    The reason d'etre for cutting objects is their overloaded subsref,
%    which allows a syntax like X.cutting(t0,t1).
%
%    See also tsig/cut.

U = CollectInStruct(X);
U = class(U, mfilename);

