function A=address(get, put, rm, figh);
% address - address object constructor
%    A = address(getter, putter, killer) creates an address object A for
%    storage of Dynamic objects.
%
%    address(A), where A is already an address object, returns A.
%    
%    See Dynamic, address/subsref, GUIaccess.

if nargin<1, % void
    [get, put, rm, figh] = deal([]);
elseif isa(get, 'address'),
    A = get;
    return;
elseif nargin==1, 
    error(['Wrong single input arg of class ''' class(get) '''.'])
end
A = CollectInStruct(get, put, rm, figh);
A = class(A,mfilename);
superiorto('dynamic'); % prevent overloading of get/put/rm by classes using address






