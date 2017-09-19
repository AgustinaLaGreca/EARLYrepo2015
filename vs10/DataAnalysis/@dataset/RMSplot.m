function RMSplot(D, Chan, varargin);
% dataset/RMSplot - plot RMS of recordings from a dataset
%    RMSplot(D, Chan) plots the STD of the mean across reps against the 
%    value of the varied parameter X. Chan is the recording channel.
%
%    RMSplot(D, Chan, ...) passes additional arguments to xplot.
%
%    See also dataset/rateplot.


GSP = GenericStimparams(D);
X = D.Stim.Presentation.X;
Y = D.Stim.Presentation.Y;
for icond=1:GSP.Ncond,
    [Y dt t0] = anamean(D, Chan, icond);
    Tw = -t0+burstwindow(D,icond);
    Y = getsnips(Y*1e3, dt, 0, Tw+t0); % V -> mV
    R(icond) = std(Y);
end

xplot(X.PlotVal, R, 'linestyle', '-', 'marker', 'o', varargin{:});
set(gca,'fontsize',10)
xlabel([lower(X.ParName) ' (' X.ParUnit ')'],'fontsize',12);
if isequal('Octave', X.PlotScale) || isequal('log', X.PlotScale),
    xlog125;
end
ylabel('rms (mV)','fontsize',10);
title(IDstring(D, 'full'), 'fontsize', 12, 'fontweight', 'bold', 'interpreter', 'none');



