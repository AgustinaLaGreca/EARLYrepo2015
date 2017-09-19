function T = probetubecalib(Probe, Cavity, Tube);
% probetubecalib - constructor for Probetubecalib objects
%   Probetubecalib(Probe, Cave, Tube) constructs a probetubecalib object 
%   from the triple of Transfer objects: 
%        Probe:  D/A to probe microphone transfer
%       Cavity:  D/A to cavity microphone transfer
%         Tube:  probe-to-cavity acoustic transfer constructed from Probe
%                and Cavity.
%
%   Probetubecalib(), with no input args, returns a void Probetubecalib
%   object. 
%
%   Probetubecalib(S), where S is a struct, attempts to convert S into  a
%   ProbeTubeCalib object.
%
%   Type "methods Probetubecalib" to know what can be done with and to 
%   Probetubecalib objects.

Tvoid = VoidStruct('Descr', 'Probe', 'Cavity', 'Tube', 'EditHistory');
if nargin<1, 
    T = Tvoid;
elseif nargin==1 & isstruct(Probe), % merely a conversion from struct to PTC object 
    T = union(Tvoid,Probe);
elseif nargin==1 & isa(Probe, mfilename), % tautology
    T = Probe;
    return;
else,
    Descr = description(Tube);
    EditHistory = [];
    T = CollectInStruct(Descr, Probe, Cavity, Tube, EditHistory);
    Descr = description(Tube);
    EditHistory = [];
end

T = class(T, mfilename);













