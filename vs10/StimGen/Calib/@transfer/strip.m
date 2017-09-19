function T = strip(T);
% transfer/strip - strip Transfer object of its contents.
%    T = strip(T) removes the actual data from T, but keeps the other
%    fields that identify it. This serves to save disk space when 
%    writing copies of to disk. The Userdata field of T can be used to
%    hold the information needed to restore the stripped content. See
%    transfer/unstrip for more details.
%
%    See also dataset/ds2trf, transfer/unstrip, transfer/setuserdata.

[T.Freq, T.Ztrf] = deal([]);



