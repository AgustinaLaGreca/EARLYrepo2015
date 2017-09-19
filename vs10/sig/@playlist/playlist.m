function P = playlist(varargin);
% playlist/playlist - playlist constructor
%   P = playlist(W1, W2, ...) returns an playlist object containing
%   waveforms W1, W2 .. The waveforms are numerical vectors having arbitrary
%   lengths. The order and rep counts of the waveforms can be  specified 
%   using playlist/list.
%
%   See also playlist/list, "methods playlist".

iplay = []; % will contain play list
Nrep = []; % will contain rep counts
if nargin<1,
    Waveform = {};
    P = CollectInStruct(Waveform, iplay, Nrep);
elseif isequal(1,nargin) && isstruct(varargin{1}),
    P = varargin{1}; % struct -> playlist conversion
elseif isequal(1,nargin) && isa(varargin{1},mfilename),
    P = varargin{1}; return; % playlit -> playlist "conversion" (return now to avoid class call @ EOF)
elseif ~all(cellfun(@isnumeric, varargin)),
    error('All playlist waveforms must be numerical arrays.');
elseif ~all(cellfun(@isvector, varargin)),
    error('All playlist waveforms must be numerical arrays.');
else, % waveforms okay; store them as column vectors
    Waveform = Columnize(varargin, 'cellwise');
    P = CollectInStruct(Waveform, iplay, Nrep);
end
P = class(P, mfilename);



