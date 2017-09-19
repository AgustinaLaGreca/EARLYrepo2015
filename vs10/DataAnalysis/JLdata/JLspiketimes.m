function SPT = JLspiketimes(Uidx)
% JLspiketimes - returns arrival times driven spikes of recording
%    JLspiketimes(Uidx) returns an array holding the arrival times (ms) of
%    spikes during the steady-state portion of the stimulus, referred to th
%    onset of the stimulus.

if ~isscalar(Uidx),
    error('Uidx arg must be single recording index.')
end

W = JLwaveforms(Uidx);
SPT = W.SPTraw(W.APinStim)-W.t_stimOnset;





