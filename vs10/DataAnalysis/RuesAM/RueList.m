function [L, N] = RueList;
% RueList - list of Newnames audiogram cell names (paired cells only)
%   [L, Nrep] = RueList returns the list in cell string L, Nrep is vector
%   containing #reps for each cell. Only p* names, i.e. cells that are one
%   of a pair.
%
%   See also RueceptiveField.

%Ddir = 'D:\Data\RueData\IBW\Audiograms\NewNames';
%Ddir = 'D:\Data\RueData\IBW\Including_Zero_dB';
Ddir = 'D:\Data\RueData\IBW2';
if isequal('Cel',CompuName),
    Ddir = 'C:\D_Drive\rawData\RueXcorr';
end

L = dir(fullfile(Ddir, 'p*_0_.ibw'));
L = {L.name};
L = strrep(L, '_0_.ibw', '');

N = zeros(1,numel(L));
for ii=1:numel(L),
    qq = dir(fullfile(Ddir,[L{ii} '*.ibw']));
    N(ii) = numel(qq);
end





