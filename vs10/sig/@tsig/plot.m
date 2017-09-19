function H = plot(T, varargin)
% tsig/plot - plot tsig object
%   plot(T) or T.plot plots T using a proper time axis. If T is a
%   multi-channel tsig, the different channels are plotted in different
%   colors. For complex waveforms, the imaginary part is plotted using a
%   dotted line. Plot returns an array of handles to the lines of the plot.
%
%   plot(T,..) passes any trailing arguments to Matlab's PLOT function.
%   Note: color specs are overruled by the default Matlab color order.
%
%   See tsig, PLOT, "methodshelp tsig".


% extract time and waveform vectors in cell array; each cell is channel
tc = time(T); % time in ms
zc = waveform(T); 
IL = all(haslogic(T)); % all channels are logical values
% cycle through Matlab default colors
CLR = get(gca,'ColorOrder');
Nc = size(CLR,1); % # colors in CLR

H = [];
for ii=1:nchan(T),
    t = tc{ii};
    z = zc{ii};
    x = real(z);
    y = imag(z);
    icol = 1 + rem(ii-1,Nc);
    if ii==1,
        h = plot(t, x, varargin{:}, 'color', CLR(icol,:));
    else,
        h = xplot(t, x, varargin{:}, 'color', CLR(icol,:));
    end
    if ~isreal(z),
        h2 = xplot(t, y, varargin{:}, 'color', CLR(icol,:));
        % select lines among h2 that do show a line
        hnonline = findobj(h2, 'linestyle', 'none');
        h2line = setdiff(h2,hnonline);
        set(h2line,'linestyle', ':');
        h = [h h2];
    end
    H = [H, h];
end
xlabel('Time (ms)');
if IL,
    ylim([-0.1 1.1]);
    ylabel('Logical value');
    set(gca,'ytick',[0 1], 'yticklabel',{'no' 'yes'})
else,
    ylabel('Amplitude (A.U.)');
end
figure(gcf);

if nargout<1, clear H; end % only return handles if requested

