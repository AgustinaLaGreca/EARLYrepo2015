function PL = JLpowerlaw(Uidx);
% JLpowerlaw - instantaneous firing probability from averaged subthr potential
%    JLpowerlaw(Uidx);

% merely book-keeping
AD = JLanovaDetails(Uidx);
CS = JLcycleStats(Uidx);

ad = structpart(AD, {'UniqueRecordingIndex' 'mean2d' 'spikeprob___________________' ...
    'NspikeStim' 'mean_at_spt' 'var_at_spt'    'Vinst'    'histV'    ... 
    'hist_mean_at_spt'    'hist_var_at_spt'    'cond_spike_prob'});
cs = structpart(CS, {'mean___________________' 'MeanSpont1' ...
    'MeanOnset'  'MeanDriv' 'MeanTail' 'MeanOffset' 'MeanSpont2'});

tim.AnaDur = CS.AnaDur;
tim.NsamPerCycle = AD.NsamPerCycle;
tim.DurPerCell = 1e-3*CS.AnaDur/AD.NsamPerCycle.^2;  % s duration represented by one cell in 2D grid

PL = structJoin(ad,cs, '-time', tim);



