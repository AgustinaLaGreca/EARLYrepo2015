function Ph= ztrf(EC, ichan, Freq);
% earcalib/ztrf - complex transfer factor of earcalib object
%    ztrf(EC, ichan, Freq) returns the complex factor of the
%          AD(ichan)-->eardrum 
%    transfer function stored in Earcalib object EC, evaluated at frequency 
%    Freq [Hz].
%
%    If waveform having a complex spectral value Cspec at Freq Hz is sent
%    to channel ichan D/A converter, the corresponding spectral component 
%    at the eardrum will be
%
%         Cspec_EAR = Cspec_DA * ztrf(EC,Freq)
%
%    This can be used to realize a requested spectrum at the aerdrum.
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
Ph = interp1(frequency(Trf), ztrf(Trf), Freq(:));
Ph = reshape(Ph, Sz);

