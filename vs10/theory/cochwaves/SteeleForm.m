function sf = SteeleForm(kL, kb);
%  SteeleForm form factor of 3D fluid motion
%     SteeleForm(kL, kb);
%     INPUTS
%       kL: k*L, with k the wavenumber and L the scala radius in compatible
%           units (e.g. rad/mm and mm, respectively)
%       kb: k*b, with k the wavenumber and b the effective BM "radius",
%       i.e., pi times the effective width of the BM.
%
%     The effective fluid mass equals SteeleForm(k*L,k*b)./k when using the
%     BM velocity as generalized velocity. It equals SteeleForm(k*L,k*b).*k
%     when using the longitudinal volume flow as generalized velocity. The
%     latter option is convenient in a transmission-line treatment of the
%     traveling wave.
%
%     See also coch3Dk.

[kL, kb] = SameSize(kL, kb);

K1L = besselk(1,kL);
K1b = besselk(1,kb);
I1L = besseli(1,kL);
I1b = besseli(1,kb);
K0b = besselk(0,kb);
I0b = besseli(0,kb);

sf = -(K1L.*I0b + I1L.*K0b)./(K1L.*I1b - I1L.*K1b);



