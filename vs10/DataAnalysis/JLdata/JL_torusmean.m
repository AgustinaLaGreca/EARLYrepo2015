function [ph1, ph2] = JL_torusmean(H, P);
% JL_torusmean(N2d, P)

[N1, N2] = size(H);
H = (H(:)-mean(H(:)));
H = H/std(H(:));
H = H.^P;
H = H/sum(H);
[Phi2, Phi1] = meshgrid((0:N2-1)/N2, (0:N1-1)/N1);
Phasor1 = exp(2*pi*i*Phi1(:));
Phasor2 = exp(2*pi*i*Phi2(:));

m1 = sum(H.*Phasor1);
m2 = sum(H.*Phasor2);

ph1 = mod(cangle(m1),1);
ph2 = mod(cangle(m2),1);




