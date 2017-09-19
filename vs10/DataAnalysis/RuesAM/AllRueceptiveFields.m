function AllRueceptiveFields;
% AllRueceptiveFields - compute all rueceptive fields & save in XLS files
%
%   See also RueceptiveField.

%Sdir = 'D:\Data\RueData\IBW\Audiograms\NewNames\Overlap\autoXLS';
Sdir = 'D:\Data\RueData\IBW\Including_Zero_dB\FRAs';

L = RueList;
for ii=1:numel(L),
    fn= L{ii}
    CC = RueceptiveField(fn);
    xlswrite(fullfile(Sdir,fn), CC);
end





