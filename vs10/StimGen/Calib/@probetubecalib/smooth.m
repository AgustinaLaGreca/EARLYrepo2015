function P = Smooth(P, Freq, BW);
% Probetubecalib/Smooth - Smooth part of Tube response 
%   PTC = Smooth(PTC, [f1 f2], BW) smooths the complex trf of P.Tube
%   over a limited frequency range of the Tube trf of PTC. BW is the given
%   bandwidth is octaves.
%
%   See also Probetubecalib/edit, Probetubecalib/clip, smoothen.

if ischar(BW), BW = eval(BW); end

ihit = find(betwixt(frequency(P.Tube), Freq(1), Freq(2)));
if isempty(ihit), return; end

Z = ztrf(P.Tube);
limZ = Z(ihit); % Z limited to requested f
offset = min(ihit)-1; % sample offset of Zlim freq within Z

df = df_oct(P.Tube); % freq spacing in octaves
[limZ, ism] = smoothen(limZ, BW, df);
% smoothly insert the smoothed part back into Z
ism = min(ism,numel(limZ));
Z = insert(Z,offset,limZ,ism);

% now put this back in P.Tube
tub = struct(P.Tube); % ugly, but for now faster than writing a transfer/smooth
tub.Ztrf = Z;
P.Tube = transfer(tub);

P.EditHistory = [P.EditHistory struct('fcn', fhandle(mfilename), 'args', {{probetubecalib, Freq, BW}})];








