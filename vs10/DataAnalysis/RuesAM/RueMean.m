function RueMean;
% RueMean - compute across-rep averages of Rue audiogram data & save

% Ddir = 'D:\Data\RueData\IBW\Audiograms\NewNames';
% Mdir = 'D:\Data\RueData\IBW\Audiograms\NewNames\means';
% Ddir = 'D:\Data\RueData\IBW\Audiograms\NewNames';
% Mdir = 'D:\Data\RueData\IBW\Audiograms\NewNames\means';
% Ddir = 'D:\Data\RueData\IBW\Including_Zero_dB';
% Mdir = 'D:\Data\RueData\IBW\Including_Zero_dB\means';
Ddir = 'D:\Data\RueData\IBW2';
Mdir = 'D:\Data\RueData\IBW2\means';
qq = dir(fullfile(Ddir,'p*_0_.ibw'));
FNS = strrep({qq.name}, '_0_.ibw', '');
for ifile=1:numel(FNS);
    local_mean(Ddir, Mdir, FNS{ifile});
end


%=============================
function local_mean(Ddir, Mdir, FN);
FN
qq = dir(fullfile(Ddir, [FN '_*.ibw']));
if isempty(qq),
    error(['Data files ''' FN ''' not found.']);
end
Nrep = numel(qq);
FNS = {qq.name};
Y = 0;
for irep = 1:Nrep, 
    fn = FNS{irep};
    D = IBWread(fullfile(Ddir, fn));
    Y = Y + D.y;
end
Y = Y/Nrep;
OutFile = fullfile(Mdir, [FN '_mean.mat']);
save(OutFile, 'Y');


