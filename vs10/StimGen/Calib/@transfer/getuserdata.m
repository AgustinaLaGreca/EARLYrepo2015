function ud = getuserdata(T, fn);
% transfer/getuserdata - get userdata of transfer object
%    getuserdata(T) returns the userdata struct of transfer object T.
%
%    getuserdata(T, 'Foo') returns filed Foo of the userdata.
%    
%    See also transfer/setuserdata.

ud = T.Userdata;
if nargin>1,
    ud = ud.(fn);
end







