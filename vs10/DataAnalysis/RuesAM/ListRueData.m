function [F, Nrep, Ncell] = ListRueData(dtype);
%  ListRueData - compile list of Rue datafiles of given type
%   [FN, Nrep, Ncell] = ListRueData or [FN, Nrep] = ListRueData('audiograms'); 
%   [FN, Nrep, Ncell] = ListRueData('FM'); 
%
%   FN is cellstring of filenames
%   Nrep(k) is # reps of cell k whose name is FN{k}.
%   Ncell is # cells.

if nargin<1, dtype = 'audiograms'; end
doPool = ismember('*', dtype) || ismember('q', dtype);
dtype = strrep(dtype, '*', '');
dtype = strrep(dtype, 'q', '');

switch upper(dtype(1)),
    case 'A',
        Ddir = 'D:\Data\RueData\IBW\Audiograms'
    case 'F',
        Ddir = 'D:\Data\RueData\IBW\FM'
    otherwise,
        error('dtype must be Audiogram or FM.');
end

qq = dir(fullfile(Ddir, ['*.ibw']));
Nfile = numel(qq);
isIpsi = false(1,Nfile);
for ifile=1:Nfile,
    fn = qq(ifile).name;
    isIpsi(ifile) = ismember('i', strtok(fn, '.'));
    FN{ifile} = strtok(fn, '_');
end

if doPool,
    FN = FN(isIpsi);
    FN = strrep(FN, 'i', 'q');
end

[F, I, J]  = unique(FN);
Ncell = numel(F);
for ii=1:Ncell,
    Nrep(ii) = sum(J==ii); % # files starting with F{ii}
end






