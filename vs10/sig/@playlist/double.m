function X = double(P,k)
% playlist/double - expand playlist into numerical vector.
%
%   Double(P) or P.double expands the waveforms in playlist P into a single
%   numerical column vector, counting repetition.
%
%   Note: use Double(~P) to ignore reps.
%
%   See playlist/not, playlist, "methods playlist".

X = [];
for ii=1:length(P.iplay),
    iwave = P.iplay(ii);
    Nrep = P.Nrep(ii);
    X = [X; repmat(P.Waveform{iwave}, Nrep, 1)];
end






