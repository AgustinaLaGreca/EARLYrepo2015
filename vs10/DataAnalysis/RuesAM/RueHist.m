function Q = RueHist(FN);
% RueHist - histograms of correlations of individual 80-dB stim conditions
%    Q = RueHist('p101022') displays histograms of


global HackHist; HackHist = 1;

CFN = mfilename;
CP = FN;

Q = getcache(CFN, CP);

FNA = [FN 'Ac'];
FNB = [FN 'Bc'];
if isempty(Q),
    global HackHist; HackHist = 1;
    RueCrossShuf(FNA, FNA);
    public kloeng;
    Q.XAA = kloeng(:,:,1:29); % 80-dB only
    global HackHist; HackHist = 1;
    RueCrossShuf(FNB, FNB);
    public kloeng;
    Q.XBB = kloeng(:,:,1:29); % 80-dB only
    global HackHist; HackHist = 1;
    RueCrossShuf(FNA, FNB);
    public kloeng;
    Q.XAB = kloeng(:,:,1:29); % 80-dB only
    putcache(CFN, 200, CP, Q);
end

local_hist(Q.XAA, FNA,'');
local_hist(Q.XBB, FNB,'');
local_hist(Q.XAB, FNA,FNB);

function local_hist(X, FN1, FN2);
fh = figure;
set(fh,'units', 'normalized', 'position', [0.117 0.148 0.86 0.755])
Freq = logispace(1,48.5,29); % freq in kHz
for ifreq=1:29,
    BinC = linspace(-1,1,20);
    ha = subplot(5,6,1+ifreq);
    x = X(:,:,ifreq);
    idiag = 1:(1+size(x,1)):numel(x);  
    Pnd = x; Pnd(idiag)=[];
    Pd = x(idiag);
    hist(Pnd, BinC);
    hh = findobj(ha,'type', 'patch');
    hold on;
    hist(Pd, BinC);
    set(hh, 'facecol', 'r');
    xlim([-1.05 1.05]);
    title([num2str(Freq(ifreq),3) ' kHz']);
    %title([P.FN1 ' / ' P.FN2]);
end
subplot(5,6,1);
text(0.1, 0.5, strvcat(FN1,FN2));
set(fh, 'PaperOrientation', 'landscape', 'PaperUnits', 'normalized', 'PaperPosition', [0.1 0.1 0.85 0.85]);
print(fh, '-djpeg', fullfile('D:\Data\RueData\IBW2\microhist\', ['microhist_' FN1 '_' FN2]));






