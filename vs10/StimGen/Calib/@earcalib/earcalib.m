function EC = earcalib(X);
% earcalib - constructor for earcalib objects
%   Earcalib() constructs a yet void Earcalib object, to be filled later by
%   Earclib/GUI. Earcalib objects contain calibration data measured in situ. 
%   
%   Earcalib(S), where S is a struct, attempts to convert S into  an 
%   Earcalib object.
%
%   Type "methods earcalib" to know what can be done with and to Earcalib 
%   objects.

ECvoid = VoidStruct('Experiment', 'iCalib', 'activeDAchan', 'Param', 'Transfer', 'EditHistory');
if nargin<1, 
    EC = ECvoid;
elseif nargin==1 & isstruct(X), % merely a conversion from struct to earcalib object 
    EC = union(ECvoid,X);
elseif nargin==1 & isa(X, mfilename), % tautology
    EC = X;
    return;
else,
    error('Wrong type of input argument. Consult Earcalib help text.');
end

EC = class(EC, mfilename);













