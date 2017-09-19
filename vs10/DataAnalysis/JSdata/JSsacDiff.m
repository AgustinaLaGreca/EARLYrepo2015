function JSplotPeaks(ExpID, iRec, warp1, warp2);
% JSplotPeak1 - plot peak changes (i.e., N+T minus N) against S/N ratio.
%   EXAMPLES
%    JSplotPeak1('RG11295',56:65) peak changes of warp=(0,0) correlograms
%    JSplotPeak1('RG11295',56:65, 0.3) peak changes of warp=(0.3,0.3) correlograms
%    JSplotPeak1('RG11295',56:65, -0.1, 0.2) peak changes of warp=(-0.1,0.2) correlograms

if nargin<3,
    warp1 = 0;
end
if nargin<4,
    warp2 = warp1; % see help
end
[MC, Stim] = JSmcoutput(ExpID, iRec);
MC = cellify(MC);
ds = read(dataset, ExpID, iRec(end));
icell = ds.ID.iCell;
allNoiseSPL = unique([Stim.MnoiseSPL]);
figure;
set(gcf,'units', 'normalized', 'position', [0.254 0.169 0.272 0.761], 'PaperPositionMode', 'auto')
NnoiseSPL = numel(allNoiseSPL);
ha = [];
for ii=1:NnoiseSPL, % select all recs w given noise SPL
    noiseSPL = allNoiseSPL(ii);
    ihit = find([Stim.MnoiseSPL]==noiseSPL);
    [mc, stim] = sortAccord(MC(ihit), Stim(ihit), -[Stim(ihit).ToneSPL]); % decreasing S/N ratio
    N_SN = numel(mc); % # S/N ratio values for this noise SPL
    sacDiff = local_xfield(mc, warp1, warp2, 'totDiff');
    [sacDiff, dum, SNR] = denan(sacDiff, [stim.ToneSPL]);
    ha(end+1) = subplot(NnoiseSPL, 1, NnoiseSPL+1-ii);
    CLR = [0.75 0 0];
    plot(SNR, sacDiff, 'v-', 'color', CLR, 'markerfacecolor', CLR, 'linewidth', 3, 'markersize', 8);
    if ii==1, 
        xlabel('S/N ratio (dB)');
    end
    ylabel('correlogram diff');
    if ii==NnoiseSPL,
        title([name(ds.Experiment) ' cell ' num2str(icell), ...
            '  warp=(' num2str(warp1) ',' num2str(warp1) ')']);
    end
    legend({['noise: ' num2str(noiseSPL) ' dB/Hz']}, 'location', 'northwest');
    legend(gca,'boxoff');
end
set(gcf,'units', 'normalized', 'position', [0.254 0.169 0.272 0.1+NnoiseSPL/4*0.661], 'PaperPositionMode', 'auto')
JSlinkaxes(ha, 'x')
%===================================
function FV = local_xfield(mc, warp1, warp2, FN);
for ii=1:numel(mc),
    mci = mc{ii};
    ihit = find(abs([mci.warp1]-warp1)<1e-3 & abs([mci.warp2]-warp2)<1e-3);
    if isempty(ihit),
        warning('Requested warp values not found in dataset');
        FV(ii) = nan;
    else,
        FV(ii) = mci(ihit).(FN);
    end
end







