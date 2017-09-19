function CC = RueceptiveField(FN, doPlot);
% RueceptiveField - receptive fields of Rue audiograms based on autocorr
%   CC = RueceptiveField(FN, doPlot)

if nargin<2, doPlot=(nargout<1); end
%Ddir = 'D:\Data\RueData\IBW\Audiograms\NewNames'; 
%Ddir = 'D:\Data\RueData\IBW\Including_Zero_dB';
Ddir = 'D:\Data\RueData\IBW2';

SPL = 0:10:80;
NSPL = numel(SPL); % 8 SPLs in audiogram
Freq = logispace(1,48.5,29); % freq in kHz
Nfreq = numel(Freq); % 29 freqs in audiogram
Ncond = Nfreq*NSPL; 

if nargin<1, % all files
    local_all(Ddir, NSPL, Nfreq);
    return;
end

% cache
CFN = mfilename;
CFP = lower(FN);
CC = getcache(CFN, CFP);
if isempty(CC), % compute it
    disp('... computing ...');
    [CC, NSPL, Nfreq] = local_computeCC(FN, Ddir);
    putcache(CFN, 1e3, CFP, CC);
end
if doPlot,
    Rueceplot(CC,-1, 1,1,1); % 1,1,1 = draw xlabels, ylabels, colorbars
    title(FN, 'fontweight', 'bold');
end


%=========================
function [CC, NSPL, Nfreq] = local_computeCC(FN, Ddir);
SPL = 0:10:80;
NSPL = numel(SPL); % 8 SPLs in audiogram
Freq = logispace(1,48.5,29); % freq in kHz
Nfreq = numel(Freq); % 29 freqs in audiogram
Ncond = Nfreq*NSPL; 
if nargin<1,
    CC = [];
    return;
end
qq = dir(fullfile(Ddir, [FN '_*.ibw']));
if isempty(qq),
    error(['Data files ''' FN ''' not found.']);
end
Nrep = numel(qq);
FNS = {qq.name};
Y = [];
for irep = 1:min(Nrep,15), % max 15 reps to prevent memory overload
    fn = FNS{irep};
    D = IBWread(fullfile(Ddir, fn));
    Y = [Y, D.y];
end
NsamRep = size(Y,1); % # samples in one rep of the whole audiogram
NsamCond = round(NsamRep/Ncond); % # samples in one stim condition

CC = zeros(Ncond,1); % vector holding suffled auto correlations per stim cond
for icond=1:Ncond,
    irange = (icond-1)*NsamCond + (1:NsamCond); % samples indices of this stim cond
    CC(icond) = shufcorr(Y(irange,:)); % one cond, all reps
end
CC = reshape(CC, [Nfreq, NSPL]).';


function local_all(Ddir, NSPL, Nfreq);
qq = dir(fullfile(Ddir, '*_0_.ibw'));
FNS = {qq.name};
Nfile = FNS;
figure; colormap(flipud(gray)); set(gcf, 'PaperPosition', [0.5 0.5 20 28]);
ScreenPos = [0.263 0.0234 0.723 0.909];
set(gcf,'units', 'normalized', 'position', ScreenPos);
isub = 1; Nsub = 24;
Nrow = 6; Ncol = 4;
Ncell = numel(FNS);
for icell = 1:Ncell,
    subplot(Nrow,Ncol, isub);
    fn = strrep(FNS{icell}, '_0_.ibw', '');
    CC = RueceptiveField(fn);
    Leftmost = isequal(1,rem(isub,Ncol));  %leftmost row of plots: provide y-tickmarks, ylabels
    Bottom = (isub>(Nrow-1)*Ncol) ...
        || (icell>Ncell-Ncol); %bottow row of plots: provide x-tickmarks, xlabels
    Rueceplot(CC, -1, Bottom, Leftmost,0); % last 0: no colorbar
    title(fn, 'fontweight', 'bold');
    if isub==Ncol*Nrow,
        figure;
        set(gcf, 'PaperPosition', [0.5 0.5 20 28]);
        ScreenPos = ScreenPos + [0.025    -0.02  0  0];
        set(gcf,'units', 'normalized', 'position', ScreenPos);
        colormap(flipud(gray));
        isub = 1;
    else,
        isub = isub + 1;
    end
    pause(0.05);
end


