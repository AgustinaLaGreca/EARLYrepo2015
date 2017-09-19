function APstats(ExpID, RecID, icond, SlopeDur, ThrSlope);
% APstats - evaluate AP statistics of a recording

if nargin<2, RecID = []; end
if nargin<3, icond = []; end
if nargin<4, SlopeDur = 0.25; end % ms default duration over which to evaluate slopes
if nargin<5, ThrSlope = -3; end % V/s default upper boundary of slopes considered

[D, DS, E] = readTKABF(ExpID, RecID, icond);

ichan = 1; % channel 1 is neural waveform
Vrec = D.AD(ichan).samples; 
dt = D.AD(ichan).dt_ms;
dVdt = diff(Vrec)/dt;
Nsam = numel(dVdt);

% compute min & max values of dVdt over a moving window of length SlopeDur
NsamSlope = round(SlopeDur/dt);
minSlope = -runmax(-dVdt, NsamSlope);

minSlope(minSlope>ThrSlope) = ThrSlope;
%dplot(dt, minSlope)
set(figure,'units', 'normalized', 'position', [0.733 0.376 0.231 0.517]);
hist(minSlope(minSlope<ThrSlope), 20);

qsteep = ismember(dVdt, minSlope);

% uminSlope = unique();
% B = linspace(min(minSlope), ThrSlope, 5);







