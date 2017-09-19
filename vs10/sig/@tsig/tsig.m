function T = tsig(Fsam, Waveform, t0, varargin);
% tsig - tsig constructor
%   tsig with no input arguments returns a void tsig object.
%
%   T=tsig(Fsam, W, t0) creates a tsig object T with sample frequency
%   Fsam  [kHz], waveform(s) W, and onset t0 [ms]. Default onset is 0 ms.
%   The waveform W is a numerical or logical vector. 
%
%   Note that t0 is a "sample offset": t0 is rounded toward the nearest 
%   integer multiple of 1/Fsam. Use fft techniques to realize subsample
%   time shifts.
%
%   tsig(Fsam, {V W ..}, t0) returns a multi-channel tsig object. When t0
%   is a vector, its elements indicate the onsets of the corresponding
%   channels.
%
%   Type T.help for tsig/subsref help; T.help = 0 for tsig/subsasgn help.
%
%   See also tsig/isvoid, "methodshelp tsig".

if nargin<1, % void tsig
    Fsam = [];
end
if nargin<2,
    Waveform = {};
end
if nargin<3,
    t0 = 0;
end
Internal = struct(varargin{:});
%-----check input args-------

% specials cases: struct input or object input
% mfilename
% class(Fsam)
if nargin==1 && isstruct(Fsam),
    T = Fsam;
    T = class(T, mfilename);
    return;
elseif nargin==1 && isa(Fsam, mfilename),
    T = Fsam;
    return;
elseif nargin<1,
    t0 = [];
    T = CollectInStruct(Fsam, Waveform, t0, Internal);
    T = class(T, mfilename);
    return;
end

%------regular call from here: tsig is fully specified------------
if nargin<2,
    error('No waveforms specified.');
end

% sample freq
FSerr = 'Fsam must be positive scalar.';
if ~isempty(Waveform) && isempty(Fsam),
    error(FSerr);
end
if nargin>0,
    if ~isscalar(Fsam),
        error(FSerr);
    end
    error(numericTest(Fsam, 'real/positive/noninf', 'Fsam is '));
end

% make sure waveform is in column vector format, & packed in cell array
if isnumeric(Waveform) || islogical(Waveform), Waveform = {Waveform}; end

if ~iscell(Waveform) ...
        || ~all(cellfun(@isnumeric, Waveform) | cellfun(@islogical, Waveform)) ...
        || ~all(cellfun(@isvector, Waveform) | cellfun(@isempty, Waveform)), % vector or [] is okay
    error('Waveform must be numerical vector or cell array of numerical vectors.');
end
Waveform = Columnize(Waveform, 'cellwise'); % reshape all elements to be col vectors
Waveform = Waveform(:).'; % row cell array containing channels in numerical col vectors 

% t0 (= onset time)
TZerr = 't0 must be real scalar.';
if ~isempty(Fsam) && isempty(t0),
    error(TZerr);
end
if ~isempty(t0),
    error(numericTest(t0, 'real/noninf/nonnan', 't0 is '))
    if ~isscalar(t0),
        t0 = t0(:).'; % row vector
    end
    try,
        [t0, Waveform] = SameSize(t0, Waveform);
    catch,
        error('Sizes of onset and Waveform do not match.');
    end
    DT = 1/Fsam; % sample period in ms
    t0 = DT*round(t0/DT); % t0 is sample offset; see help text
end


T = CollectInStruct(Fsam, Waveform, t0, Internal);
T = class(T, mfilename);

