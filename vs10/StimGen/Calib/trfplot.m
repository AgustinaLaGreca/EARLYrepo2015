function trfplot(C, plotArg);
% trfplot - quick & dirty transfer function plotter
%   trfplot(C, plotArg);
%   trfplot({C1 C2}, plotArg);

if nargin<2, plotArg='n'; end

if iscell(C), % transfer C1 to C2
    C{1}.TRF = C{2}.TRF./C{1}.TRF;
    comp_lag = C{1}.Lag_ms - C{2}.Lag_ms;
    comp_phase = 1e-3*C{1}.Freq.*comp_lag;
    phasor = exp(2*pi*i*comp_phase);
    C = C{1};
    C.TRF = C.TRF.*phasor;
end

subplot(2,1,1); 
xplot(C.Freq/1e3, A2dB(abs(C.TRF)), plotArg);
ylabel('Magnitude (dB)')
grid on

subplot(2,1,2); 
xplot(C.Freq/1e3, cunwrap(angle(C.TRF)/2/pi), plotArg);
xlabel('Frequency (kHz)')
ylabel('Phase (cycles)');
grid on


