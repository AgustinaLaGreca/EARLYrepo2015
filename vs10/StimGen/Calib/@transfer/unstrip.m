function T = unstrip(T, flag);
% transfer/unstrip - restore previously stripped contents of transfer object
%    T = unstrip(T) attempts to restore content previously removed by
%    transfer/strip. In order for this to work, the userdata field of T has
%    to provide the information needed to rebuild T the same way t was
%    originally constructed. It must have a field 'recreate' holding a cell
%    array C such that T = feval(C{:}) recomputes T.
%    See dataset/ds2trf for an example. 
%
%    T = unstrip(T, 'force') also recomputes T if T already holds data. The
%    default behavior is to return T unchanged if it already holds data.
%
%    See also dataset/ds2trf, transfer/unstrip, transfer/setuserdata.

[flag] = arginDefaults('flag', '');

if isempty(T.Ztrf) || isequal('force', lower(flag)), % stripped or force
    UD = getuserdata(T);
    T = feval(UD.recreate{:});
end



