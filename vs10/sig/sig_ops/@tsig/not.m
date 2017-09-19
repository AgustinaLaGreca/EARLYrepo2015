function S = not(S);
% tsig/not - NOT for tsig objects
%    not(S) performs a sample-wise, channel-wise logical NOT on
%    tsig object S.
%
%    See also tsig/plus, tsig/relop.

S.Waveform = cellfun(@not, S.Waveform, 'UniformOutput',false);



