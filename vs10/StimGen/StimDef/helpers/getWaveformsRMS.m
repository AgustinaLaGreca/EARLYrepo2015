function Y= getWaveformsRMS(W)
% Waveform/samples - extract rms value from serie of waveforms object
%     [Y dt] = RMS(W) 
%
%     [Y dt] = samples(W, 'all') compute the rms value of the signal
%     including silence chunk.

[Nwf, Nchan] = size(W);
Y = zeros(Nwf,Nchan);
for iwave=1:Nwf
    Y(iwave,:) = W(iwave,:).RMS;
end
    
end

