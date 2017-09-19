function Lstr = CM_IO(igerbil, irec, Delay, plotArg);
% CM_IO = I/O characteristics of single-tone CM responses
%   CM_IO(igerbil, irec, Delay, plotArg);
%
%   EXAMPLE
%   CM_IO(78,1+[1297 1595 1279 1399], [40 540 210 -7]);

if nargin<3, Delay=0; end % 0 us default Delay to be compensated
if nargin<4, plotArg=''; end % default: new plot
if isempty(plotArg), initplot = @cla; else, initplot=@nope; end

if numel(irec)>1 || numel(Delay)>1, % multiple recordings: recursion
    [irec, Delay] = SameSize(irec, Delay);
    for ii=1:numel(irec),
        Lstr{ii} = CM_IO(igerbil, irec(ii), Delay(ii), ploco(ii));
    end
    f1;
    subplot(2,3,1);
    legend(Lstr{:}, 'location', 'northeast');
    return;
end

%========single recording from here============
S = ZWOAEimport(igerbil, irec);
if ~isequal(0, S.Fsingle) || ~isequal(0, S.Fsup) || ~isequal(1, S.Nzwuis),
    error('Stimulus is not a single tone.');
end
if ~isequal('CM',S.RecType),
    warning('Not a CM recording.');
end

% extract all harmonics of stim freq below Nyquist
Fnyq = S.fs/2; % Nyquist freq in kHz
Nmax = floor(Fnyq/S.Fzwuis); % # harmonics considered
Freq = (1:Nmax).'*S.Fzwuis;
% restrict spectrum to Freq components. Add zero DC term.
MG = [0; SpecComp(S.df, S.MG, Freq)];
PH = [0; SpecComp(S.df, S.PH, Freq)];
PH = cunwrap(PH);
Nspec = 2^10; 
Cspec = dB2A(MG).*exp(2*pi*i*PH);
Ncomp = numel(Cspec);
Tsig = real(ifft([Cspec; zeros(Nspec-Ncomp,1)]));
dph = 1/Nspec;
% extract stimulus phase from acoustical recording
SpAc = ZWOAEsubspec(igerbil, S.companionID, 'zwuis');
PHref_rad = 2*pi*(SpAc.PH + S.PHzwuis);
%======plot===========
f1;
set(gcf,'units', 'normalized', 'position', [0.17563 0.44167 0.80438 0.48417]);
% Magn
subplot(2,3,1); initplot();
xplot([0;Freq],MG,['.-' plotArg])
xlim([0 Fnyq]);
ylim(max(MG)+[-85 5]);
Lstr = [num2str(S.Lzwuis) ' dB SPL'];
legend(Lstr, 'location', 'northeast');
ylabel('Magnitude (dB)');
% Phase
subplot(2,3,4); initplot();
xplot([0;Freq],PH,['o-' plotArg])
xlabel('Frequency (kHz)');
ylabel('Phase (Cycle)')
% Waveform
dphase = PHref_rad-2*pi*Delay*1e-3*S.Fzwuis; % delay in us, Fzwuis in kHz -> phase shift in rad
Stim = dB2A(S.Lzwuis)*cos(dphase+2*pi*(0:Nspec-1)/Nspec).';
subplot(2,3,[2 6]); initplot();
xdplot(dph, [Tsig; Tsig], ['-' plotArg]);
xdplot(dph, max(abs(Tsig))*[Stim; Stim]/max(abs(Stim)), [':' plotArg]);
xlabel('Phase (Cycle)');
f2;
[dum, imax] = max(Stim); idesc = 1+rem(imax+round(0.25*Nspec)-1,Nspec);
SC = 1;%db2a(-S.Lzwuis);
xplot(SC*Stim,SC*Tsig,plotArg);
xplot(SC*Stim(idesc),SC*Tsig(idesc),['o' plotArg]);
xplot(SC*Stim(idesc+5),SC*Tsig(idesc+5),['*' plotArg]);
xlabel('Instantaneous stim Ampl');
ylabel('Instantaneous resp Ampl');




