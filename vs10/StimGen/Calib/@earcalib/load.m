function [EC, FFN] = load(EC, XP, ifile);
% Earcalib/load - load Earcalib object from file.
%    EC = load(Earcalib(), Exp, K) reads a the k-th Ear calibration file of
%    experiment Exp into an Earcalib object EC. Default K is inf, which by 
%    convention means the most recent one. Default Exp is the current 
%    Experiment. Note that the first input is a dummy serving to Methodize
%    this function.
%
%    load(Earcalib(), 'Foo', ...) uses find(experiment, 'Foo') as
%    Experiment.
%
%    load(Earcalib(), Exp, NaN) prompts the user to select one of the
%    calibrataion files of the experiment.
%
%    [EC, FFN] = load(...) also returns the full filename where EC was
%    retrieved.
%
%    See also Methodizing, Earcalib/save, Experiment.

if nargin<2,
    XP = [];
end
if nargin<3,
    ifile = inf; % most recent
end
if ischar(XP),
    XP = find(experiment, XP);
end

if isequal([],XP),
    XP = current(experiment);
end

if isvoid(XP),
    error('Experiment not found or void.');
end

FileFilter = fullfile(folder(XP), [name(XP) '_*' fileExtension(EC)]);
if isnan(ifile),
    [FN, PP]=uigetfile(FileFilter, 'Select Ear calibration file.');
    if isequal(0,FN), FFN=''; return; end;
    FFN = fullfile(PP, FN);
elseif ifile<inf,
    FFN = strrep(FileFilter, '*', dec2base(ifile,10,3));
else, % most recent file
    FF = dir(FileFilter);
    FF = FF(~[FF.isdir]);
    if isempty(FF),
        error(['No ear calibration data found for experiment ''' name(XP)  '''.']);
    end
    [dum, irecent] = max([FF.datenum]);
    nam = FF(irecent).name;
    FFN = fullfile(folder(XP), nam);
end
if ~exist(FFN,'file'),
    error(['Ear calib file ''' FFN  ''' not found.']);
end

load(FFN,'EC','-mat');





