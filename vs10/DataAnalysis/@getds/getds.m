function G = getds(XP);
% getds - construct Getds object
%   DS = getds('RG11306') returns a "dataset grabber" for experiment RG11306.
%   When used with an index, DS will return a dataset:
%   DS(2) 
%   ans =
%   RG11306, Rec 2  <1-2-FS> Pen 0, NaN um  Rotterdam/CLUST 20-Jun-2011 16:18:08
%      FS stim   11 x 20 x 10/50 ms     (complete)
%      Carrier frequency = 10000:20000 Hz (log steps)
%         RW_electrode:  10999 ms @ 111.6 kHz (1227607 samples)
%
%   See getds/subref for constructing pooled_dataset objects using curly
%   braces.

[XP, Nmax] = arginDefaults('XP/Nmax', [],100); % Nmax = max # cached datasets

G = CollectInStruct(XP, Nmax);
G = class(G, mfilename);




