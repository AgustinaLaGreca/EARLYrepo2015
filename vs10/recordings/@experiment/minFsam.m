function Fm = minFsam(E, Dev);
% experiment/minFsam - minimum sample rate of experiment
%   minFsam(E) returns the minimum sample rate in Hz dictated by the 
%   specifications of E entered during Experiment/edit.
%   recordings specified for experiment E. minFsam returns -inf if no
%   analog recordings are specified for E. This includes the case of a void
%   experiment.
%
%   minFsam(E, Dev), where Dev is one of 'RX6', 'RP2', returns minumum
%   sample rate of the named recording device (default: RX6).
%
%   See also MaxFsam, Experiment.

if nargin<2, Dev = 'RX6'; end

if isvoid(E), Fm = -inf; return; end; 

Fm.RX6 = -inf;
Fm.RP2 = -inf;

% Audio
dev = E.Audio.Device;
Fm.(dev) = max(Fm.(E.Audio.Device), E.Audio.MinFsam_kHz);
% Recording
if canrecord(E),
    dev = E.Recording.General.Device;
    Fm.(dev) = max(Fm.(dev), E.Recording.General.MinFsam);
end
% Electra
if canelectrify(E),
    dev = E.Electra.Device;
    Fm.(dev) = max(Fm.(dev), E.Electra.MinFsam_kHz);
end

Fm = Fm.(Dev);
fcn = fhandle([Dev 'sampleRate']); % RX6samplerate, etc
Fall = fcn('all');
Fm = min(Fall(Fall>=Fm));








