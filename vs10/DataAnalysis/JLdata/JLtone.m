function S = JLtone(ExpID, RecID, icond, AP, Cdelay, plotArg);
% JLtone - periodicty analysis of tone sweeps from JL's recordings
%    JLtone(ExpID, RecID, icond, AP, plotArg);

[AP, Cdelay, plotArg] = arginDefaults('AP/Cdelay/plotArg', [], 7.5, 'n');

AP = replaceMatch(AP, [], -10); % default -10 V/s AP thr slope
if numel(AP)<3,
    AP = [AP(1) -1 1.4]; %truncation 0.35 ms before peak to 1.00 ms after peak
end
if numel(AP)<4, % cut-out limit equal to acceptance limit
    AP = [AP AP(1)];
end
ThrAPdownRate = AP([1 4]); % first thr is acceptance as AP; second (smaller abs) thr for being cut out
APwindow = AP(2:3); % AP truncation window

if iscell(RecID),
    for ii=numel(RecID):-1:1,
        S(ii) = JLtone(ExpID, RecID{ii}, icond, AP, Cdelay, ploco(ii));
        legend({'Both' 'Contra' 'Ipsi' });
        ylim auto;
        DS = sgsrdataset(ExpID, RecID{end});
        Freq = round(DS.xval(icond));
        if numel(RecID)>2,
            title([ExpID '  ', RecID{1} '/' RecID{2} '/' RecID{3} '  ' num2str(Freq(1)) ' Hz  ' num2str(max(DS.SPL)) ' dB SPL' ], 'fontweight', 'bold');
        end
    end
    fenceplot(S(ii).SPTraw, ylim, 'color', 0.7*[1 1 1])
    return;
elseif numel(icond)>1,
    '*****'
    for ii=1:numel(icond),
        S(ii) = JLtone(ExpID, RecID, icond(ii), AP, Cdelay, ploco(ii));
    end
    return;
end

MyFlag('LastRecordingAnalyzed', {ExpID, RecID, icond});
MaxOrd = 7;

AvMethod = 'mean'; % default
% if ischar(Vclip),
%     [AvMethod, Mess] = keywordmatch(lower(Vclip), {'mean', 'median'}, 'flag');
%     error(Mess);
% end

Ttrans1 = 13; % ms first part during stim ignored
Ttrans2 = 0;  % ms last part during stim ignored

[D, DS, L]=readTKABF(ExpID, RecID, icond);
if ~isequal('FSLOG', DS.stimtype) && ~isequal('FS', DS.stimtype),
    error([DS.title ' is not a FSLOG dataset.']);
end

TTT = [num2str(icond) ': ' num2str(DS.xval(icond)) ' ' DS.x.Unit '  ' DS.info];

SPL = DS.SPL;
ichan = find(SPL>=0);

Fsam = 0.999993*D.AD(1).Fsam_kHz; % kHz (% correct for different clock rates of A/D and D/A)
dt = 1/Fsam; % sample period in ms
burstDur = DS.burstdur(1); % ms tone duration
postDur = DS.repdur-burstDur; % ms post-tone dur
Freq1 = DS.fcar(icond,1); Period1 = 1e3/Freq1; % cycle in ms
Freq2 = DS.fcar(icond,2); Period2 = 1e3/Freq2; % cycle in ms
if ~isequal(Freq1,Freq2),
    error('Unequal left & right frequencies.');
end
Nsam1 = round(1e3*Fsam/Freq1);
Tcar = 1e3/Freq1;
Nsamdelay = round(1e3*Fsam/Cdelay);
% stim1 = D.AD(3).samples; % ipsi stimulus waveform
% stim2 = D.AD(4).samples; % contra stimulus waveform
rec = cat(2, D.AD(2:end,1).samples); % recording; ignore first sweep (silent)
Nrep = size(rec,2);
[SPTraw, Tsnip, Snip, APslope] = deal([]);
for irep=1:Nrep,
    [rec(:,irep), tSPTraw, Tsnip, tSnip, tAPslope] = APtruncate2(rec(:,irep), ThrAPdownRate, dt, APwindow);
    SPTraw = [SPTraw; tSPTraw];
    Snip = [Snip, tSnip];
    APslope = [APslope tAPslope];
end
AC = shufcorr(rec); % shuffled autocorr across reps
stdrec = std(rec,[], 2); % std across reps
rec = mean(rec,2); % average over reps
% truncate spikes

tmid = 0.5*(Ttrans1 + burstDur-Ttrans2); % midddle of burst
Trec = burstDur-(Ttrans1+Ttrans2); % max duration of useful portion
NcarCycle = floor(Trec/Tcar) % # carrier cycles
NsamCycle = round(Tcar/dt) % # samples per carrier cycle
imid = 1 + round(tmid/dt); % index of middle sample of analysis window
NsamRec = NsamCycle*NcarCycle % total # cycles in analysis window
isam = imid - round(NsamRec/2) + (1:NsamRec); % sample range of analysis window
medRec = median(rec);

% average across cylces
[rec, stdrec] = deal(rec(isam), stdrec(isam));
rec = rec-mean(rec);
%NsamCycle, NsamRec
Y = LoopMean(rec, NsamCycle);
stdY = LoopMean(stdrec, NsamCycle);
ioffset = isam(1)-1 + Nsamdelay;
Y = circshift(Y,ioffset);
stdY = circshift(stdY,ioffset);
Y = [Y;Y]; % 2 cycles
stdY = [stdY;stdY]; % 2 cycles

S =CollectInStruct(ExpID, RecID, icond, '-', ...
    Freq1, Freq2, SPL, '-', ...
    Cdelay, ThrAPdownRate, APwindow, Ttrans1, Ttrans2, '-', ...
    dt, '-', SPTraw, Tsnip, Snip, APslope, ...
    '-', Y, stdY, AC ...
    );

% plot waveforms
xdplot(dt, Y, plotArg);
set(gcf,'units', 'normalized', 'position', [0.0188 0.462 0.98 0.438]);
%ylim(ylim);
xlabel('Time (ms)');
TracePlotInterface(gcf);
Thalf = dt*NsamCycle;
title(TTT)


















