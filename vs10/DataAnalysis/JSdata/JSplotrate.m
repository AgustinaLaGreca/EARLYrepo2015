function JSplotRate(ExpID, iRec, warp);
% JSplotRate - plot spike rate increase (N+T versus N) against S/N ratio.
%   EXAMPLES
%    JSplotRate('RG11295',56:65) rate plot(s) of warp=0 recordings
%    JSplotRate('RG11295',56:65, 0.3) rate plot(s) of warp=0.3 recordings

if nargin<3,
    warp = 0;
end

[MC, Stim, CMX] = JSmcoutput(ExpID, iRec);
MC = cellify(MC);
ds = read(dataset, ExpID, iRec(end));
icell = ds.ID.iCell;
allNoiseSPL = unique([Stim.MnoiseSPL]);
figure;
set(gcf,'units', 'normalized', 'position', [0.512 0.172 0.272 0.761], 'PaperPositionMode', 'auto')
NnoiseSPL = numel(allNoiseSPL);
ha = [];
for ii=1:NnoiseSPL, % select all recs w given noise SPL
    noiseSPL = allNoiseSPL(ii);
    ihit = find([Stim.MnoiseSPL]==noiseSPL);
    [cmx, stim] = sortAccord(CMX(ihit), Stim(ihit), -[Stim(ihit).ToneSPL]); % decreasing S/N ratio
    N_SN = numel(cmx); % # S/N ratio values for this noise SPL
    meanRateN = local_xfield(cmx, warp, 'meanRateN');
    meanRateT = local_xfield(cmx, warp, 'meanRateT');
    stdRateN = local_xfield(cmx, warp, 'stdRateN');
    stdRateT = local_xfield(cmx, warp, 'stdRateT');
    RateDiff = stdRateT - stdRateN;
    stdRateDiff = sqrt(stdRateT.^2 + stdRateN.^2);
    [RateDiff, dum, SNR] = denan(RateDiff, [stim.ToneSPL]);
    ha(end+1) = subplot(NnoiseSPL, 1, NnoiseSPL+1-ii);
    CLR = [0 0.5 0];
    plot(SNR, RateDiff, 's-', 'color', CLR, 'markerfacecolor', CLR, 'linewidth', 3, 'markersize', 8);
    hold on;
    errorbar(SNR, RateDiff, stdRateDiff, 'color', CLR, 'linewidth', 2)
    if ii==1, 
        xlabel('S/N ratio (dB)');
    end
    ylabel('Rate increase (sp/s)');
    if ii==NnoiseSPL,
        title([name(ds.Experiment) ' cell ' num2str(icell), ...
            '  warp=' num2str(warp)]);
    end
    legend({['noise: ' num2str(noiseSPL) ' dB/Hz']}, 'location', 'northwest');
    legend(gca, 'boxoff')
end
set(gcf,'units', 'normalized', 'position', [0.512 0.172 0.272 0.1+NnoiseSPL/4*0.661], 'PaperPositionMode', 'auto')
JSlinkaxes(ha, 'x')
%===================================
function FV = local_xfield(mc, warp, FN);
for ii=1:numel(mc),
    mci = mc(ii);
    ihit = find(abs([mci.warp]-warp)<1e-3);
    if isempty(ihit),
        warning('Requested warp values not found in dataset');
        FV(ii) = nan;
    else,
        FV(ii) = mci.(FN)(ihit);
    end
end







