function D = Zweig(omega,x,delta)
%  Zweig - stapes to BM transfer according to Zweig et al, 1976
%     D = Zweig(omega,x,delta)
%       x: normalized distance from stapes, 0<x<1

N = 5; % total # wave cycles to cutoff point
Frat = 1000; % ratio between smallest and largest resonance freq
[omega,x,delta] = SameSize(omega,x,delta);
om_r = Frat.^(1-x)/sqrt(Frat); % resonance freq. Equals 1 halfway the cochlea by convention
max(om_r), min(om_r)
D = i.*omega.*om_r.^(-3/2).*exp(-4*i*N*asin(omega./om_r-i*delta/2))./(om_r.^2-(omega-i*delta.*om_r/2).^2).^(3/4);





