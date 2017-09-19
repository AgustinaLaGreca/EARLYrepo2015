function Ph= Phase(EC, ichan, Freq);
% earcalib/magnitude - magnitude of earcalib object
%    magnitude(EC, ichan, Freq) returns the magnitude [dB] of the
%          AD(ichan)-->eardrum 
%    transfer function stored in Earcalib object EC, evaluated at frequency 
%    Freq [Hz].
%
%    If a tone of Freq Hz is sent to the D/A converter having a magnitude
%    of magn_DA dB, the level at the eardrum will be
%
%         SPL_EAR = magn_DA + magnitude(EC,Freq)  dB
%
%    This can be used to realize a requested SPL at the eardrum.
%    Freq may be an array. All its elements must lie in the frequency range
%    of EC.
%
%    See Earcalib, Earcalib/GUI, Earcalib/frequency.

if ~isfilled(EC),
    error('Transfer object EC is not filled.');
end

if any(~infreq_range(EC,Freq,ichan)),
    error('Freq exceeds frequency range of Earcalib.')
end
Trf = EC.Transfer(ichan);

Sz = size(Freq);
Ph = interp1(frequency(Trf), magnitude(Trf), Freq(:));
Ph = reshape(Ph, Sz);

