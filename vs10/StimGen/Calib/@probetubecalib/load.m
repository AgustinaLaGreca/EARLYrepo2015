function [T, FFN] = load(T, FN, Prompt);
% Probetubecalib/load - load Probetubecalib object from file.
%    T = load(probetubecalib(),'Foo') reads a transfer object from file
%    Foo.transfer in the calibration dir that was previosuly saved using
%    Transfer/save. Explicit directory or extension will overrule the 
%    default dir datadir(PTC) and extension fileExtension(PTC), where
%    PTC is some Probetubecalib object.
%
%    Note that the first argument is a dummy for Methodizing this function.
%
%    T = load(T0,'?',P) prompts for the file name using UIGETFILE
%    with prompt P. The default prompt is 'Specify file name of
%    Probetubecalib object'. If the user cancels, T0 is returned.
%
%    T = load(T0,'?Foo', ...) also prompts the user, but this time using a
%    filter Foo.
%
%    [T, FFN] = load(...) also returns the full file name FFN where the
%    data were retrieved.
%
%    See also Methodizing, Probetubecalib/save, UIGETFILE.

if nargin<2,
    FN = '?';
end
if nargin<3,
    Prompt = 'Specify file name of Probetubecalib data.';
end
FFN = '';

if isequal('?', FN(1)), % prompt for file
    FileFilter = FN(2:end);
    if isempty(FileFilter),
        FileFilter = '*';
    end
    FileFilter = FullFilename(FileFilter, datadir(T), fileExtension(T));
    [FN, PATH]  = uigetfile(FileFilter, Prompt);
    if isequal(0, FN), return; end % user cancelled
    FN = fullfile(PATH, FN);
end

FFN = FullFilename(FN, datadir(T), fileExtension(T));
warning off MATLAB:elementsNowStruc
load(FFN,'T','-mat');
warning on MATLAB:elementsNowStruc
T = local_fix(T);
[Dum fn  dum] = fileparts(FFN);
T.Descr = [fn ': ' T.Descr];


%===============
function T = local_fix(T);
% repair backward compatibility of transfer objects within T
if isstruct(T.Probe),
    T.Probe = transfer(T.Probe);
end
if isstruct(T.Cavity),
    T.Cavity = transfer(T.Cavity);
end
if isstruct(T.Tube),
    T.Tube = transfer(T.Tube);
end






