function D = datadir(dum, D);
% Probetubecalib/datadir - get or set directory for Probetubecalib data
%   datadir(Probetubecalib()) returns the local directory for storing
%   ProbeTubeCalib data. If none is defined, an error occurs.
%   
%   datadir(Probetubecalib(), 'd:\Foo\subdir') sets the ProbeTuubeCalib
%   directory to d:\FOO\subdir. The specified directory must exist. This 
%   setting is rembered across Matlab sessions.
%
%   NOTE 
%   The first input arg is a dummy for Methodizing this function.
%
%   See also Datadir, setuplist, Methodizing.

if nargin<2, % get
    D = FromSetupFile('localsettings', 'Calibrationdir', '-default', 123);
    if isequal(123,D), 
        error('No probetubecalib directory defined. Use probetubecalib/datadir to specify it.');
    end
else,  % set
    if ~isdir(D),
        error(['''' D ''' is not an existing directory.'])
    end
    ToSetupFile('localsettings', '-propval', 'Calibrationdir', D);
end




