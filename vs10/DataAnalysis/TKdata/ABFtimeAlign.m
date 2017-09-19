function D = ABFtimeAlign(D);
% ABFtimeAlign - correct misalignment across reps ("sweeps") of ABF recording
%   D = ABFtimeAlign(D) corrects D for timing misaligment due to
%   incommensurate sample rates of AD and DA.
%
%   It is assumed that the recording is in chan 1, that the stimulus is 
%   in chan 2, and that the first 20 ms of the stimulus contains enough
%   non-zero samples to be used for alignment. Alignment resolutions
%   is limited to the sample period of the ABF data (usually 10 us).
%
%   See also readTFABF, TKpool.

ichanRec = 1; 
if D.Nchan==2,
    ichanStim = 2; 
elseif D.Nchan==4,
    ichanStim = [3 4];
else,
    error(['Don''t know how to time align datasets having ' num2str(D.Nchan) ' recording channels.']);
end

dt = D.dt_ms; % sample period in ms

if D.Nsweep<3, return; end % single stim sweep: nothing to align
% collect first 20 ms of stim in matrix; each col is sweep.
NsamStim = round(50/dt); % # samples ~ 50 ms
NsamMaxShift = round(2/dt); % max 2 ms shift

for ichan=ichanStim(:).',
    stim = [];
    for isweep = 2:D.Nsweep,
        sam = D.AD(isweep,ichan).samples(1:NsamStim);
        stim = [stim, sam];
    end
    Avstim(:,ichan) = mean(stim,2); % mean stim acros sweeps
end
% select stim channel having strongest signal
[dum, ichan] = max(var(Avstim));
Avstim = Avstim(:,ichan);
%dplot(dt, Avstim); pause
HWin = hann(NsamStim); % apply hann window to suppress boundary effects when crosscorrelating
Avstim = Avstim.*HWin;

for isweep = 2:D.Nsweep,
    [MaxC, iLag(isweep-1)] = maxcorr(HWin.*stim(:,isweep-1), Avstim, NsamMaxShift);
end
iLag = iLag-min(iLag);
NsamNew = size(D.AD(1,1).samples, 1) - max(iLag);

iLag = [0 iLag]; % add entry for 1st sweep

% apply the alignment
irange = 1:NsamNew;
for isweep = 1:D.Nsweep,
    samShift = iLag(isweep);
    for ichan = 1:D.Nchan,
        D.AD(isweep,ichan).samples ...
            = D.AD(isweep,ichan).samples(samShift+irange);
    end
end
D.NsamPerSweep = NsamNew;
D.sweepDur_ms = D.NsamPerSweep*dt;






