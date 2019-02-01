function [DACstr, Lstr, Rstr] = getDACstr(AudioChannel, RecordingSide)
% DACSTR - Returns DACstr to be used on the GUI for the DAC field
%   The values depend on the audio channel that is being used and the
%   recording side (experiment field).

%   See Experiment, getrecordingSideStr

% Figure out the recording side
[Lstr, Rstr] = getRecordingSideStr(RecordingSide);

switch AudioChannel
    case 'Left', DACstr = {Lstr};
    case 'Right', DACstr = {Rstr};
    case 'Both', DACstr = {Lstr Rstr 'Both'};
end
end