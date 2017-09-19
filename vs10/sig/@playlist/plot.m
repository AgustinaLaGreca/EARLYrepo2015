function h = plot(P, varargin);
% playlist/PLOT - PLOT playlist object.
%
%   PLOT(P) plots playlist object P. A list must be specified for P.
%   A line graph is produced displaying the expanded waveform of P, reps
%   included. The start of each next waveform is indicated by a plot symbol.
%   PLOT returns a 2-element array containing the handle to the line and
%   symbols, respectively.
%
%   PLOT(P,...) passes any remaining args to PLOT, for instance color, 
%   linestyle, etc.
%
%   PLOT(X,P,...) plots the samples of P against numerical vector X.
%
%   Note: use PLOT(~P) to display reps only once.
%
%   See playlist/double, "methodshelp playlist".

if isnumeric(P), % % must be plot(X,P, ..) syntax; see help
    if ~isPlaylist(varargin{1}),
        error('Either first or second argument of tsig/plot must be tsig object.');
    end
    X = P;
    P = varargin{1};
    varargin = varargin(2:end);
else, % no X axis specified
    X = 1:nsamplay(P);
end

if isempty(P.iplay),
    error('No list specified for P. Use playlist/list.')
end
Y = double(P);
h = plot(X, Y, varargin{:});
c = get(h,'color');
istart = Offsets(P); % position of starts
h = [h xplot(X(istart), Y(istart), 'o', 'color', c)];
figure(gcf);
if nargout<1, clear h; end % only return handle when requested




