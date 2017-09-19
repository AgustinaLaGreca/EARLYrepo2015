function A=address(D, Ax, Nam);
% Dataset/address - get/set address of uploaded dynamic object
%   A=address(D) returns the address of dynamic object D, i.e. a struct A
%   with function_handle-valued fields get, set, rm that provide access to
%   D as follows:
%         X=A.get() returns the uploaded value of D;
%         A.put(D2) changes the uploaded value to D2;
%         A.rm() removes the uploaded version of D.
%
%   D=address(D,A) sets the adress to A. This has the side effects of
%   uploading D to its new address and removing any previous storage of D.
%
%   D=address(D,figh,'Foo') is the same as
%   D=address(D,GUIaccess(figh,'Foo'))
%
%   See Dynamic/upload, dynamic/hasbeenuploaded, GUIaccess.

if nargin<2, % get
A = D.access;
else, % set
    if nargin>2, Ax = GUIaccess(Ax,Nam); end
    % remove any earlier storage
    if ~isempty(D.access), D.access.rm(); end
    % store new address and upload D to it
    D.access = Ax;
    D = upload(D);
    A = D; % return new D
end

