function Ph= Phase(EC, ichan, Freq, Tau);
% earcalib/phase - phase of earcalib object
%    phase(EC, ichan, Freq) returns the phase [Cycle] of the
%          AD(ichan)-->eardrum 
%    transfer function stored in Earcalib object EC, evaluated at frequency 
%    Freq [Hz]. Importantly, this phase does *not* include the phase lag
%    originating from the "wideband delay" WBdelay(EC,ichan).
%
%    If a tone of Freq Hz is sent to the D/A converter having a phase of 
%    phase_DA cycles, the phase at the eardrum will be
%
%         phase_EAR = phase_DA + Phase(EC,Freq)
%
%    This can be used to realize a requested phase at the aerdrum.
%    Freq may be an array. All its elements must lie in the frequency range
%    of EC.
%
%    phase(EC, ichan, Freq, Tau) uses a wideband delay of Tau ms instead of
%    EC's own wideband delay. For the special values Tau='flat', and Tau =
%    'wflat', see Transfer/setWBdelay.
%
%    phase(EC, ichan, Freq, Tau) uses a wideband delay of Tau ms instead of
%    EC's own wideband delay.
%
%    See Earcalib, Earcalib/GUI, Earcalib/frequency, Earcalib/WBdelay.

if nargin<4,
    Tau = [];
end

if ~isfilled(EC),
    error('Transfer object EC is not filled.');
end

if any(~infreq_range(EC,Freq,ichan)),
    error('Freq exceeds frequency range of Earcalib.')
end
Trf = EC.Transfer(ichan);
if isempty(Tau),
    deltaTau = 0;
else,
    deltaTau = Tau-getWBdelay(Trf);
end
if ~isequal(0,deltaTau),
    Trf = setWBdelay(Trf,Tau);
end

Sz = size(Freq);
Ph = interp1(frequency(Trf), phase(Trf), Freq(:));
Ph = reshape(Ph, Sz);










