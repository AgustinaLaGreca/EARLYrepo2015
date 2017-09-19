function A=afteryou(A,varargin);
% action/after - 	 interdependent ready state
%    afteryou(A, Ax1, Ax2, ..) specifies a list of action objects A1, A2, ...
%    with pointers Ax1, Ax2, .. whose current status determines the being 
%    ready of A itself. isready(A) will return True if all of the listed
%    objects (A1,A2..) are in the 'finished' state. This is a kind of 
%    "dynamic overloading". See action/isready for details.
%
%    Note that stopping ('interrupting') is not the same as finishing.
%    Stopped objects are not 'finished', and if any of the listed A1, A2,...
%    is stopped, isready will be wrong. In this case A will have to be
%    actively stopped, too. This is typically handled by GIOacion, which
%    vists all the action object of a certain GUI.
%
%    afteryou(A, figh, 'Foo', 'Faa' ..) is the same as
%    afteryou(A, GUIaccess(figh, 'Foo'), GUIaccess(figh, 'Faa') ...)
%
%    Defining an AfterYou list for A is an alternative to providing an
%    explicit subclass(A)/isready. Unless the overloaded method
%    subclass(A)/isready explicitly uses the AfterYou property of A (for
%    instance by invoking isready(action(A))), the AfterYou list will be
%    ignored if a subclass(A)/isready method exists.
%
%    See also action/isready, GUIaction, action/start, action/wrapup.

[A, Mess] = criticalDownload(A, 'initialized', 'define ''isready'' dependencies of');
error(Mess);

figh = [];
for ii=1:nargin-1,
    arg = varargin{ii};
    if isnumeric(arg),
        figh=arg;
    elseif ischar(arg), % arg is name; last figh must be GUI handle
        A.AfterYou = [A.AfterYou, GUIaccess(figh,arg)];
    elseif isstruct(arg), % address
        A.AfterYou = [A.AfterYou, arg];
    else,
        error('Wrong type of input arguments. See help text.');
    end
end
upload(A);






