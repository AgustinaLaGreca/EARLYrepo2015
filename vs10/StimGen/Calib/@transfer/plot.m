function [lh, h] = plot(T, varargin);
% transfer/plot - plot Transfer object.
%    h = plot(T) plots the amplitude and phase of the transfer function in
%    Transfer object T, and returns handles to the two lines.
%
%    plot([ah1 ah2], T) uses axes handles ah1 and ah2 to place the graphs of
%    the magnitude and phase, respectively.
%
%    plot(T,...) and plot([ah1 ah2],T, ...) pass additional arguments to
%    the builtin PLOT command used to realize the graphs.

if ~isa(T, 'transfer'), % plot(ah, T, ...) syntax. Rearrange input args
    if nargin<2 || ~isa(varargin{1}, 'transfer'),
        error('Either the first or the second input arg of Transfer/Plot must be Transfer object.');
    end
    if any(~istypedhandle(T,'axes')) || ~isequal(2,numel(T)),
        error('First input arg of Transfer/Plot must be either a Transfer object or a pair of axes handles.');
    end
    % now rearrange the input args
    [ah1, ah2] = DealElements(T);
    T = varargin{1};
    varargin = varargin(2:end);
else, % plot(T, ...) syntax. Provide own axes handles.
    set(gcf,'units', 'normalized', 'position', [0.324 0.436 0.637 0.478]);
    ah1 = subplot(2,1,1);
    ah2 = subplot(2,1,2);
    set(ah1, 'NextPlot', 'add');
    set(ah2, 'NextPlot', 'add');
end
h = [ah1 ah2];

if ~isfilled(T),
    error('Transfer object T is not filled.');
end


axes(ah1);
lh(1)=plot(T.Freq/1e3, A2dB(abs(T.Ztrf)), varargin{:}); 
grid on
xlog125; 
ylabel('Amp (dB)');
title(T.Description);
hold on;

axes(ah2);
lh(2)=plot(T.Freq/1e3, cunwrap(cangle(T.Ztrf)), varargin{:}); 
grid on
ylabel('Phase (cycle)'); 
xlabel('Frequency (kHz)');
title(['Wideband delay (not shown): ' num2str(T.WB_delay_ms) ' ms']);
hold on;




