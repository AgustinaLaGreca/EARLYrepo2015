function [CC, meanCC, isPair, CCall, DD] = RueCrossShuf(FN1, FN2, doPlot);
% RueCross - shuffled condition-wise crosscorr of means of audiograms
%   [CC, SS, isPair] = RueCross(FN1, FN2);
%
%   Example:
%   RueCross('p90206Ac', 'p90206Bc')

if nargin<3, doPlot=(nargout<1); end

CFN = mfilename; CacheParam = {FN1, FN2};
CC = getcache(CFN, CacheParam);
if isempty(CC), % try swapping args
    CC = getcache(CFN, {FN2, FN1});
end
global HackHist;
if isempty(HackHist), HackHist= 0; end
if isempty(CC) || HackHist, % compute it
    SPL = 0:10:80;
    NSPL = numel(SPL); % 8 SPLs in audiogram
    Freq = logispace(1,48.5,29); % freq in kHz
    Nfreq = numel(Freq); % 29 freqs in audiogram
    Ncond = Nfreq*NSPL;

    [RL, Nrep] = RueList;
    i1 = strmatch(upper(FN1),upper(RL));
    if isempty(i1), error([FN1 ' not found.']); end
    i2 = strmatch(upper(FN2),upper(RL));
    if isempty(i2), error([FN2 ' not found.']); end

    Nrep1  =Nrep(i1);
    Nrep2  =Nrep(i2);
    if HackHist, kloeng = nan(Nrep1,Nrep2, Ncond); end
    CC = []; % Npair x Ncond correlation matrix
    DD = [];
    CO1 = 0;
    for irep1=1:Nrep1,
        D1 = local_read_rep(FN1, irep1);
        %CO1 = CO1 + local_COV(D1.y,D1.y,Ncond);
        for irep2=1:Nrep2,
            D2 = local_read_rep(FN2, irep2);
            %COV2(irep2,1:Ncond) = local_COV(D2.y,D2.y,Ncond);
            cc = local_CC(D1.y,D2.y,Ncond);
            if HackHist, for icond=1:Ncond, kloeng(irep1,irep2, icond) = cc(icond); end; end;
            if irep1==irep2, 
                DD = [DD; cc];
            else,
                CC = [CC; cc];
            end % skip simultaneous pairs
        end
    end
    CCall = CC;
    % average correlation over pairs
    CC = mean(CC) ;
    % reshape CC into Freq x SPL grid
    CC = reshape(CC,Nfreq,NSPL).';
    putcache(CFN, 5e3, CacheParam, CC);
end
if HackHist, Public(kloeng); end;
clear global HackHist
Sdir = 'D:\Data\RueData\IBW2\FCCA'; 
if isequal('Cel',CompuName),
    Sdir = 'C:\D_Drive\rawData\RueXcorr\FCCA';
end
SFN = [FN1, '_' FN2 '.xls'];
xlswrite(fullfile(Sdir,SFN), CC);
meanCC = mean(CC(:));
isPair = isRuePair(FN1, FN2);

if doPlot,
    Rueceplot(CC, -1);
    title([FN1 '/' FN2], 'fontweight', 'bold');
end

%===============================================
function D = local_read_rep(FN, irep);
% Ddir = 'D:\Data\RueData\IBW\Including_Zero_dB';
Ddir = 'D:\Data\RueData\IBW2';
if isequal('Cel',CompuName),
    Ddir = 'C:\D_Drive\rawData\RueXcorr';
end
FFN = fullfile(Ddir, [FN '_' num2str(irep-1) '_.ibw']);
D = IBWread(FFN);

function cc = local_CC(y1, y2, Ncond);
NsamCond = round(numel(y1)/Ncond);
ir = 1:NsamCond;
for icond = 1:Ncond,
    offset = (icond-1)*NsamCond;
    cc(1,icond) = corr(y1(offset+ir), y2(offset+ir));
end

function cc = local_COV(y1, y2, Ncond);
% non-normalized version
NsamCond = round(numel(y1)/Ncond);
ir = 1:NsamCond;
for icond = 1:Ncond,
    offset = (icond-1)*NsamCond;
    qq = cov(y1(offset+ir), y2(offset+ir));
    cc(1,icond) = qq(2,1);
end





