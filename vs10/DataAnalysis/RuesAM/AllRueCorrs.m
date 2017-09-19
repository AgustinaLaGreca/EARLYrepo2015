function AllRueCorrs(Ncond, TimeWin);
% AllRueCorrs - crosscorr analysis among Rue's three test sets
%  

if nargin<1, Ncond=29*8; end

Ddir = 'D:\Data\RueData';

FileName = fullfile(Ddir,['RueCor' num2str(Ncond) '_' num2str(TimeWin)]);


RueCor.AA = RueCorr('A', 'A', Ncond, TimeWin);
RueCor.BB = RueCorr('B', 'B', Ncond, TimeWin);
RueCor.ZZ = RueCorr('Z', 'Z', Ncond, TimeWin);
RueCor.AB = RueCorr('A', 'B', Ncond, TimeWin);
RueCor.BZ = RueCorr('B', 'Z', Ncond, TimeWin);
RueCor.ZA = RueCorr('Z', 'A', Ncond, TimeWin);

save(FileName, 'RueCor');















