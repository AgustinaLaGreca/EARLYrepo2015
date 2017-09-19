function Y = JLspontAna(S, Fh, doPlot);
% JLspontAna - xcorr & spectral analysis spontaneous activity preceding stimulation of JLbeat data
%    Y = JLspontAna(S, Fh, doPlot) 
%   Y has fields: ID, Tlag, XC, df, Pspec.
%    
%    
%    See also JLspont

magLag = 10; % ms max lag of xcorr
doFilter = 1; % reject 622-Hz band

if nargin<2, Fh=[]; end
if nargin<3, doPlot = (nargout<1); end

if isempty(Fh), Fh = 20; end % 20 Hz default

CFN = mfilename;
CFP = {[S.UniqueRecordingIndex], Fh, magLag};
Y = getcache(CFN, CFP);
if isempty(Y),
    Y = local_compute(S, Fh, magLag, doFilter);
    putcache(CFN, 1e3, CFP, Y);
end

if doPlot,
    TT = [S(1).ExpID ' cell ' num2str(S(1).icell) '  ('  num2str(numel(S))  ' x 500 ms)']
    fh = figure;
    set(fh,'units', 'normalized', 'position', [0.13 0.506 0.795 0.396]);
    subplot(1,2,1);
    title(TT);
    plot(Y.Tlag, Y.XC);
    xplot(Y.Tlag, mean(Y.XC,2), 'k', 'linewidth', 4);
    xlabel('Lag (ms)');
    ylabel('Autocorrelation coeff');
    title(TT);
    subplot(1,2,2);
    freq = timeaxis(Y.Pspec, Y.df);
    plot(freq, Y.Pspec);
    [dum, i622] = min(abs(freq-0.622));
    [dum, i834] = min(abs(freq-0.834));
    i50 = 1+round((1:20)*0.050/Y.df);
    xplot(freq(i50), Y.Pspec(i50), 'gd');
    xplot(freq(i622), Y.Pspec(i622), 'rs' ,'markerfacecolor', 'r');
    xplot(freq(i834), Y.Pspec(i834), 'mv' ,'markerfacecolor', 'm');
    title(TT);
    xlog125([0.06 5]);
    xlabel('Frequency (kHz)');
    ylabel('Magnitude (dB)');
end


%================
function [Y, TT] = local_compute(S, Fh, magLag, doFilter);
% the real work
[SP, dt] = JLspont(S(1), Fh, doFilter);
%powerspec(dt,SP); pause;
Nmaxlag = round(magLag/dt);
Nsam = size(SP,1);
Ncond = numel(S);
df = 1/(dt*Nsam); % freq spacing in kHz
XC = [];
Pspec = 0;
for icond =1:Ncond,
    sp = JLspont(S(icond), Fh, doFilter);
    [xc, ilag] = xcorr(sp,Nmaxlag,'coeff');
    XC = [XC, xc(:)];
    Pspec = Pspec + abs(fft(sp(:).*hann(numel(sp)))).^2;
end
Pspec = P2dB(Pspec/Ncond);
Tlag = ilag*dt;
ID = structpart(S, {'ExpID'    'RecID'    'seriesID'    'iGerbil'    'icond'    'icell' '-' , ...
    'UniqueCellIndex'    'UniqueSeriesIndex'    'UniqueRecordingIndex'    'ABFname'    'ABFdir'    'JLcomment'    'MHcomment'})
Y = CollectInStruct(ID, Tlag, XC, df, Pspec);









