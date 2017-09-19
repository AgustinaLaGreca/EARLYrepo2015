function [R, alpha] = vectorstrength(spt, freq, TW);
% VectorStrength - vector strength of spike train
%   VectorStrength(Spt, Freq), where Spt is an array of spike arrival
%   times [ms], and Freq is a frequency [Hz], returns the vector strength
%   at frequency Freq of the spike times in Spt. 
%
%   VectorStrength(spt, freq, [t0 t1]) only uses spike times between t0 and
%   t1 ms. The phases are weighted inversely proportional to the number of 
%   times they were visited by the wrapping frequency. 
%
%   [R, alpha] = VectorStrength(spt, freq) returns complex R and 
%   confidence level according to Rayleigh statistics
%
%   If Spt is a cell array, R and alpha will be column arrays, with
%    [R(k) alpha(k)] = VectorStrength(spt{k});
%
%   If Spt is a cell array and Freq is an array having the same lenght,
%   then Freq(k) is used to compute the VectorStrength from Spt{k}.
%
%   If Spt is a N x M cell array and Freq is a length-N vector, then the
%   cells Spt(k,:) are pooled ("pooling across repetitions") before the
%   calculation of vector strength(s).
%
%   See also AnWin.

if nargin<3, TW = []; end
if numel(TW)==1,
    TW = [0 TW];
end

if iscell(spt), % handle by recursion (see help text)
    if isvector(spt),
        Ncond = numel(spt);
    else, % Ncond x Nrep  matrix
        Ncond = size(spt,1);
        % pool across reps
        for icond=1:Ncond,
            qq{icond} = [spt{icond,:}];
        end
        spt = qq;
    end
    % compute R & alpha for each condition
    freq = SameSize(freq(:), (1:Ncond).');
    for icond=1:Ncond,
        [R(icond,1), alpha(icond,1)] = vectorstrength(spt{icond}, freq(icond), TW);
    end
    return;
end

if ~isempty(TW), % apply analysis window
    spt = AnWin(spt,TW);
end

freq = 1e-3*freq; % Hz -> kHz, i.e. reciprocal to ms
Nspike = numel(spt);

if Nspike<2;,
   R = 0;
   alpha = 1;
   return;
end
ph = spt*freq; % phase of spike times, in cycle re period of stim

W = 1; % default: all spikes have trivial weight factor
if ~isempty(TW), % apply correction
   W = 0*ph + 1; % default weight factors of 
   ph0 = TW(1)*freq; % start phase in cycle
   ph1 = TW(2)*freq; % end phase in cycle
   Ncycle = floor(ph1-ph0); % # complete stim cycles in window
   qtail = (ph-ph0>Ncycle); % which spikes occur in tail: after last complete cycle
   W(qtail) = Ncycle/(Ncycle+1); % weight factors
   W = W/mean(W); % normalize the weight factors
end

R = mean(W.*exp(2*pi*i*ph));
if nargout>1,
   alpha = RayleighSign(R,Nspike);
end



