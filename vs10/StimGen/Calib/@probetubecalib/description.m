function Y = description(T, D);
% probetubecalib/description - get or set description of a probetubecalib object
%    description(P) returns the description of probetubecalib object P. 
%
%    P=description(P,'Foo') sets P's description to Foo.



if nargin<2, % get
    Y = T.Descr;
else, % set
    if nargout<1,
        error('Setting a description requires an output arg.');
    end
    if ~ischar(D)
        error('Description must be character string.');
    end
    T.Descr = D;
    Y = T;
end