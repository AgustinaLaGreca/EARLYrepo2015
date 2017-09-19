function S = stopguiaction(figh, Name, You, DT);
% stopguiaction - construct object that can stop all action of a GUI
%    S=stopguiaction(figh, Name, 'Foo', DT) creates an Action object that 
%    checks a stop flag every DT ms. If the flag is set, it visits all 
%    Action objects listed under handle figh and stops them. This 
%    construction is needed to avoid multiple, competing stop threads.
%    'Foo' is the name of the Action object in the same GUI (that is, 
%    having address GUIaccess(figh, 'Foo')) whose being finished implies
%    S to be finished (see action/AfterYou). Default DT is 113 ms.
%
%    Note that stopguiaction/isready does not exist. This necessitates the
%    use of AfterYou, i.e., the finishing of a stopguiaction object must be
%    made dependent on other Action object(s).
%
%    See also stopguiaction/setflag, action/AfterYou.

if nargin<4, DT = 213; end; 

S.flag = 0; % no stop requested yet
S.figh = figh; % needed to pass to GUIactionList
S = class(S, mfilename, action('initialized', DT));
S = upload(S, figh, Name);
afteryou(S, figh, You);



