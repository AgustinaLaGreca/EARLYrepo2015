function [SPT, Nspike] = AnWin(DS, W);
% Dataset/AnWin - sorted spike times within analysis window
%    SPT = AnWin(DS, [t0 t1]) returns spike times of dataset DS in a cell 
%    array SPT (See dataset/spiketimes). Only spike times falling between
%    t0 and t1 ms (re stimulus onset) are returned.
%
%    AnWin(DS, []) and AnWin(DS) are the same as AnWin(DS, [0 max(DS.burstdur)]), 
%    i.e. the maximum burstdur (note that burst duration can differ between 
%    the ears)
%
%    See also Dataset/dotraster, Dataset/AnWin.

if isvoid(DS),
    error('Void dataset.');
end

if nargin<2, W = [0 max(DS.Stim.BurstDur)]; end

SPT = spiketimes(DS);
SPT  = AnWin(SPT, W);
Nspike = cellfun('length',SPT);










