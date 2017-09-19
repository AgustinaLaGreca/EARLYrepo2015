function Fnew = tritonize(Forig, Fbase);
% tritonize - convert uberzwuis frequencies to 3-tone stim freqs
%    Fnew = tritonize(Forig, Fbase) "rounds" uberzwuis primnaries Forig
%    toward multiples of Fbase in such a way that all 3-member subsets of
%    Fnew are non-equidistant.
%    
%    Forig, Fbase, Fnew all have the same units.

Ncomp = numel(Forig);
nDF = (max(Forig)-min(Forig))/Fbase/(Ncomp-1); % normalized mean freq spacing
Dspacing = 1:Ncomp-1; % change in spacing
Dspacing = Dspacing - round(mean(Dspacing));
spacing = round(nDF) + Dspacing;
newFmin = Fbase*round(min(Forig)/Fbase);

Fnew = newFmin+Fbase*cumsum([0 spacing]);
