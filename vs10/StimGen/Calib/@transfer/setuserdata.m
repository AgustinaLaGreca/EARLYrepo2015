function T = setuserdata(T, fn, X);
% transfer/setuserdata - set userdata of transfer object
%    T = setuserdata(T, 'Foo', X) assigns teh value X to sets field Foo of
%    the userdata of transfer object T.
%
%    See also transfer/getuserdata.


if nargout<1,
    error('Too few output args. Correct syntax is: T = setuserdata(T, ''Foo'', X).');
end
T.Userdata.(fn) = X;







