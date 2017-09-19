function Delay=calibWBdelay(E, Chan);
% experiment/calibWBdelay - wideband delays of calibration data
%   calibWBdelay(E) returns the wideband delays [ms] of the calibration 
%   data of Experiment object E in a 2-element array [tau_1 tau_2]. If only 
%   one audio channel is used in E (say channel 1), the delay equals 
%   [tau_1 tau_1] by convention.
%
%   calibWBdelay(E, C) returns only the value for channel C. Channels may
%   be specified either by numbers (1, 2) or by char strings (Left, Right).
%
%   See also Transfer/getWBdelay, calibrate.

persistent PreviousValue

if nargin<2, Chan=0; end % both channels per convention

if ischar(Chan),
    [Chan, Mess] = keywordMatch(Chan, {'Left' 'Right'}, 'Audio channel');
    error(Mess);
    Chan = strmatch(Chan, {'Left' 'Right'});
end

CacheParam = {name(E) modified(E)};
% check if cached value is okay
if isstruct(PreviousValue) && isequal(CacheParam, PreviousValue.CacheParam),
    Delay = PreviousValue.Delay;
else, % get delay from data & store in cache
    Delay = localGetdelays(E);
    PreviousValue = CollectInStruct(CacheParam, Delay);
end
if ~isequal(0,Chan),
    Delay = Delay(Chan);
end

%===========================
function Delay = localGetdelays(E);
switch E.Audio.CalibrationType,
    case 'Flat',
        Delay = [0 0];
    case 'Probe',
        Probe = E.Audio.Probes;
        switch E.Audio.DAchannelsUsed,
            case 'Left',
                Delay = [1 1]*getWBdelay(cavityTrf(load(probetubecalib, Probe{1})));
            case 'Right',
                Delay = [1 1]*getWBdelay(cavityTrf(load(probetubecalib, Probe{2})));
            case 'Both',
                Delay(1,1) = getWBdelay(cavityTrf(load(probetubecalib, Probe{1})));
                Delay(1,2) = getWBdelay(cavityTrf(load(probetubecalib, Probe{2})));
        end
    otherwise,
        error(['Invalid calib type ''' E.Audio.CalibrationType '''.']);
end






