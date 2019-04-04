function Y= RMS(W, flag);
% Waveform/samples - extract rms value from waveform object
%     [Y dt] = RMS(W) returns sample array Y and sample period dt [ms]
%     from waveform. Teh waveform is fully expanded, that is, chunks are
%     repeated if specified. W may be a left/right pair, in which case Y
%     is a Nx2 matrix, but not an array.
%
%     [Y dt] = samples(W, 'all') ignores any repeats and just
%     concatenates single chunks.

[flag] = arginDefaults('flag', '');

if size(W,1)>1, error('Cannot plot multiple waveforms unless L/R pair'); end

Nchan = size(W,2);
Y = zeros(1,Nchan);
for ichan=1:Nchan
    x=[];
    for ichunk=1:numel(W(ichan).Samples)
        chunk = W(ichan).Samples{ichunk};
        if isequal('all', flag)
            x = [x; chunk];
        else
            if any(chunk)
                x = [x; chunk];
            end
        end
    end
    xrms = rms(x);
    Y(1,ichan) = xrms;
end
    
end

