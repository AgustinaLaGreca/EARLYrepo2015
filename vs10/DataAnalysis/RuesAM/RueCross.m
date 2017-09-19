function [CC, SSC, SSA] = RueCross(FN1, FN2);
% RueCross - condition-wise crosscorr of means of audiograms
%   [CC, SS] = RueCross(FN1, FN2);
%   
%   Example:
%   RueCross('p90206Ac', 'p90206Bc')

SPL = 0:10:80;
NSPL = numel(SPL); % 8 SPLs in audiogram
Freq = logispace(1,48.5,29); % freq in kHz
Nfreq = numel(Freq); % 29 freqs in audiogram
Ncond = Nfreq*NSPL; 

Y1 = local_load(FN1, Ncond);
[Y2, NsamCond] = local_load(FN2, Ncond);

CC = zeros(Ncond,1);
for icond=1:Ncond,
    CC(icond) = corr(Y1(:,icond),Y2(:,icond));
end
CC = reshape(CC,Nfreq,NSPL).';
SSC = sum(CC(:).^2)/Ncond;
AC1 = RueceptiveField(FN1);
AC2 = RueceptiveField(FN2);
SSA = sum(AC1(:).*AC2(:))/Ncond;
if nargout<1,
    Rueceplot(CC, -1);
    title([FN1 '/' FN2], 'fontweight', 'bold');
end

%===============================================

function [Y, NsamCond] = local_load(FN, Ncond);
%Mdir = 'D:\Data\RueData\IBW\Audiograms\NewNames\means';
%Mdir = 'D:\Data\RueData\IBW\Including_Zero_dB\means';
Mdir = 'D:\Data\RueData\IBW2\means';
Y = load(fullfile(Mdir, [FN '_mean.mat'])); 
Y = Y.Y;
Nsam = numel(Y);
NsamCond = round(Nsam/Ncond);
Y = reshape(Y, [NsamCond Ncond]);







