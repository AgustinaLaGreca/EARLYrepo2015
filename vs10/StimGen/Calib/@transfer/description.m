function Y = description(T, D);
% transfer/description - get or set description of a transfer object
%    description(T) returns the description of transfer object T. The
%    description of "primitive" transfer objects (filled by measurement)
%    are provided by Transfer/measure. 
%
%    T=description(T,'Foo') sets the description to Foo.


if ~isfilled(T),
    error('Transfer object T is not filled.');
end

if nargin<2, % get
    Y = T.Description;
else, % set
    if nargout<1,
        error('Setting a description requires an output arg.');
    end
    if ~ischar(D)
        error('Description must be character string.');
    end
    T.Description = D;
    Y = T;
end