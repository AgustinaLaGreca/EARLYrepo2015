function T = ds2trf(DS, Chan, Cdelay);
% dataset/ds2trf - construct transfer function from a bunch of datasets
%    ds2trf(DS, DAchan, Cdelay), where DS is an array of FS dataset or a
%    single RF or ZW dataset, constructs a Transfer object for which the stimulus 
%    tones are the input and the measured responses in DA channel DAchan are the
%    response. In overlapping regions, the dataset obtained with the
%    largest SPL takes preference. The transfer functions are derived from
%    the output of dataset/magn_phase_plot.
%
%    ds2trf(dataset(), getuserdata(T)), where T is a transfer object
%    previously computed by ds2trf, recomputes T. This is used by
%    transfer/unstrip.
%  
%    See also dataset/magn_phase_plot, transfer, transfer/unstrip.

if all(isvoid(DS)) && isstruct(Chan),
    T = local_recompute(Chan);
    return;
end

Nds = numel(DS);
if isequal('FS', stimtype(DS(1))),
    MP = Magn_Phase_Plot(DS, Chan, 'Cdelay', Cdelay);
elseif isequal('RF', stimtype(DS(1))),
    MP = [];
    for ii=1:Nds,
        mp = Magn_Phase_Plot(DS(ii), Chan, 'Cdelay', Cdelay);
        MP = [MP mp];
    end
elseif isequal('ZW', stimtype(DS(1))),
    AA = apple(DS, Chan, 0);
    MP = struct('Freq', AA.Fprim, 'Magn_dB', AA.Gain, ...
        'Phase_cycle', AA.Phase, 'iDataset', AA.irec);
else,
    error(['Datasets must have stimulus type FS, RF or ZW.']);
end
[Fmin, Fmax] = minmax(cat(1, MP.Freq)); % freq boundaries in Hz
Nfreq = round(2e3*log2(Fmax/Fmin)); % 2000 cmp per octave
Freq = logispace(Fmin, Fmax, Nfreq); % new freq axis

Z = [];; SPL = [];
for ids=1:numel(MP),
    if isequal('RF', stimtype(DS(1))),
        spl = unique(MP(ids).Yval);
    elseif isequal('FS', stimtype(DS(1))),,
        spl = DS(ids).Stim.SPL;
    elseif isequal('ZW', stimtype(DS(1))),,
        spl = DS(1).Stim.SPL;
    end
    gain = MP(ids).Magn_dB-spl+94; % normalize to Ref_resp/Pa
    ph = MP(ids).Phase_cycle;
    %ph = delayphase(ph, Freq, Cdelay/1e3, 1);
    % z = db2a(gain).*exp(2*pi*i*ph);
    GAIN = interp1(MP(ids).Freq, gain, Freq);
    PH = interp1(MP(ids).Freq, ph, Freq);
    Z = [Z; dB2A(GAIN).*exp(2*pi*i*PH)];
    SPL = [SPL; SameSize(spl, Freq)];
end
maxSPL = max(SPL,2);
hasMaxSPL = SPL==maxSPL;
Z(~hasMaxSPL) = nan;
Z = nanmean(Z,1);

% ===prepare putting stuff in transfer object
AD = anachan(DS(1), Chan); % AD data for getting specs
[dum, respunit] = conversionfactor(AD);
WB_delay = timelag(AD) + Cdelay;
% descriptors
Q_stim = 'Pressure';
Q_resp = dataType(AD);
if isequal('Microphone-1', Q_resp) || isequal('Microphone-1', Q_resp),
    Q_resp = 'Pressure';
end
Ref_stim = '1 Pa';
Ref_resp = ['1 ' respunit];
T = transfer(Q_stim, Q_resp, Ref_stim, Ref_resp);
T = ztrf(T, Freq, Z);
T = description(T, ['transfer 1 Pa -> ' Ref_resp ...
    ' constructed from ' expname(DS(1)) ' recs ' , ...
    trimspace(num2str([MP.iDataset]))]);
T = setWBdelay(T, i*WB_delay);
% Store the info need to recompute (see help text, local_recompute & transfer/unstrip)
crearg = struct('Expname', expname(DS(1)),  'iRec', irec(DS), ...
    'Chan', Chan, 'Cdelay', Cdelay);
T = setuserdata(T, 'recreate', {fhandle(mfilename) dataset() crearg});

% ========================
function T = local_recompute(UD);
DS = getds(UD.Expname);
DS = DS(UD.iRec);
T = ds2trf(DS, UD.Chan, UD.Cdelay);





