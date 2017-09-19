function T = load(T, FN, Prompt);
% transfer/load - Load Transfer object from file.
%    T = load(transfer(),'Foo') reads a transfer object from file
%    Foo.transfer in the calibration dir that was previosuly saved using
%    Transfer/save. Explicit directory or extension will overrule the 
%    default dir (calibrationdir) and extension ('transfer').
%    Note that the first argument is a dummy Transfer object, which just
%    serves to make load a Transfer method.
%
%    T = load(T0,'?',P) prompts for the file name using UIGETFILE
%    with prompt P. The default prompt is 'Specify file name of Transfer
%    object'. If the user cancels, T0 is returned.
%
%    T = load(T0,'?Foo', ...) also prompts the user, but this time using a
%    filter Foo.
%
%    See also Calibrationdir, Transfer/save, UIGETFILE.

if nargin<3,
    Prompt = 'Specify file name of Transfer data.';
end

if isequal('?', FN(1)), % prompt for file
    FileFilter = FN(2:end);
    if isempty(FileFilter),
        FileFilter = '*';
    end
    FileFilter = FullFilename(FileFilter, calibrationdir, '.transfer');
    [FN, PATH]  = uigetfile(FileFilter, Prompt);
    if isequal(0, FN), return; end % user cancelled
    FN = fullfile(PATH, FN);
end

FFN = fullfilename(FN, calibrationdir, '.transfer');

load(FFN,'T','-mat');





