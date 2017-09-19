function JSplotPeak1(ExpID, iRec, warp1, warp2);
% JSplotPeak1 - plot mean peak change (i.e., N+T minus N) against S/N ratio.
%   EXAMPLES
%   MC = JSplotPeak1('RG11295',56:65) peak changes of warp=(0,0) correlograms
%   MC = JSplotPeak1('RG11295',56:65, 0.3) peak changes of warp=(0.3,0.3) correlograms
%   MC = JSplotPeak1('RG11295',56:65, -0.1, 0.2) peak changes of warp=(-0.1,0.2) correlograms

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
NnoiseSPL = numel(allNoiseSPL);
for ii=1:NnoiseSPL, % select all recs w given noise SPL
    noiseSPL = allNoiseSPL(ii);
    ihit = find([Stim.MnoiseSPL]==noiseSPL);
    [mc, stim] = sortAccord(MC(ihit), Stim(ihit), -[Stim(ihit).ToneSPL]); % decreasing S/N ratio
    N_SN = numel(mc); % # S/N ratio values for this noise SPL
    peakDiff = local_xfield(mc, warp1, warp2, 'peak2Diff');
    [peakDiff, dum, SNR] = denan(peakDiff, [stim.ToneSPL]);
    subplot(NnoiseSPL, 1, NnoiseSPL+1-ii);
    plot(SNR, peakDiff, '*-b');
    xlabel('S/N ratio (dB)');
    ylabel('Peak2(N+T) - Peak2(N)');
    title([name(ds.Experiment) ' cell ' num2str(icell) ' ---  noise: ', ...
        num2str(noiseSPL) ' dB/Hz  warp=(' num2str(warp1) ',' num2str(warp1) ')']);
end



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







