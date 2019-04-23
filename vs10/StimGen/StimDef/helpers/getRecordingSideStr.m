function [Lstr, Rstr] = getRecordingSideStr(Recordingside)
% RECORDINGSIDESTR - Returns Lstr and Rstr as strings to be used on the GUI
%   Their values depend on the recording side that is being used, as
%   RecordingSide field in the experiment object.

%   See Experiment

switch Recordingside
    case 'Left', Lstr = 'Left=Ipsi'; Rstr = 'Right=Contra';
    case 'Right', Lstr = 'Left=Contra'; Rstr = 'Right=Ipsi';
end
end