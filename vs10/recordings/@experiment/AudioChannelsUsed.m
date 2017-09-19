function n=AudioChannelsUsed(E);
% experiment/AudioChannelsUsed - audio channels used in Experiment
%   AudioChannelsUsed(E) returns the specification of the Audio channels as
%   specified during Experiment/edit. The result is one of 'Left',
%   'Right', 'Both'.
% 
%   See also Experiment/edit.

n = E.Audio.DAchannelsUsed;
if isempty(n),
    n = '';
end






