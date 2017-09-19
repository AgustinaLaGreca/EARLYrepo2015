function [CC, meanCC] = RueCrossShuf_diagonal(FN1, doPlot);
% RueCrossShuf_diagonal - condition-wise crosscorr; version with diagonal 
%   [CC, meanCC] = RueCrossShuf_diagonal(FN1, doPlot);
%   True pairs only.
%
%   Example:
%   [CC, meanCC] = RueCrossShuf_diagonal('p90206Ac', 1)

if iscell(FN1),
    for ii=1:numel(FN1),
        [CC{ii}, meanCC(ii)] = RueCrossShuf_diagonal(FN1{ii}, doPlot);
    end
    return;
end

if nargin<3, doPlot=(nargout<1); end

FN2 = strrep(FN1, 'A' ,'B');
isPair = isRuePair(FN1, FN2);
if ~isPair,
    error('Not a true pair');
end

CFN = mfilename; CacheParam = {FN1, FN2};
CC = getcache(CFN, CacheParam);
if isempty(CC), % try swapping args
    CC = getcache(CFN, {FN2, FN1});
end
if isempty(CC), % compute it 
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
    CC = []; % Npair x Ncond correlation matrix
    for irep1=1:Nrep1,
        D1 = local_read_rep(FN1, irep1);
        for irep2=1:Nrep2,
            % following line commented out -> the diagonal IS included now
            %if irep1==irep2, continue; end % skip simultaneous pairs
            D2 = local_read_rep(FN2, irep2);
            cc = local_CC(D1.y,D2.y,Ncond);
            CC = [CC; cc];
        end
    end
    % average correlation over pairs
    CC = mean(CC) ;
    % reshape CC into Freq x SPL grid
    CC = reshape(CC,Nfreq,NSPL).';
    putcache(CFN, 5e3, CacheParam, CC);
end
Sdir = 'D:\Data\RueData\IBW\Including_Zero_dB\FCCA_diag'; 
SFN = [FN1, '_' FN2 '.xls'];
xlswrite(fullfile(Sdir,SFN), CC);
meanCC = mean(CC(:));

if doPlot,
    Rueceplot(CC, -1);
    title([FN1 '/' FN2 '  DIAG!'], 'fontweight', 'bold');
end

%===============================================
function D = local_read_rep(FN, irep);
Ddir = 'D:\Data\RueData\IBW\Including_Zero_dB';
FFN = fullfile(Ddir, [FN '_' num2str(irep-1) '_.ibw']);
D = IBWread(FFN);

function cc = local_CC(y1, y2, Ncond);
NsamCond = round(numel(y1)/Ncond);
ir = 1:NsamCond;
for icond = 1:Ncond,
    offset = (icond-1)*NsamCond;
    cc(1,icond) = corr(y1(offset+ir), y2(offset+ir));
end





