function [Spt, P, X] = peakpicker(dt, X, Thr, Reject50Hz, minISI);
% peakpicker - find peaks in recorded trace
%    T = peakpicker(dt, X, Thr, R50Hz, MinIntv)
%    Inputs:
%            dt: sample period in ms. two-component dt is [dt t0], where t0
%                is the time of the first sample (default t0=0).
%             X: recording trace. Matrices are treated as sets of rows.
%           Thr: threshold. Peaks must exceed Thr
%         R50Hz: whether to reject 50 Hz [false]
%        minISI: minimum interval (ms) between peaks (default 0.5 ms). Any
%                peak preceded or followed, within minISI ms, by a larger 
%                peak is eliminated from the output.
%   Output T is cell array. T{k} holds peak times of trace X(:,k)
%   in ms from start of X.
%
%   [T, P, Xf] = peakpicker(dt, X, ...) also returns the 
%   peak values P themselves and the 50-Hz filtered version Xf of X. P has
%   the same format as T, that is, P{k} contains the peak values extracted
%   from X(:,k). If no 50-Hz filtering was applied, Xf is identical to X.

[Reject50Hz, minISI] = arginDefaults('Reject50Hz, minISI', false, 0.5);

[X, isRow] = TempColumnize(X);
if numel(dt)>1,
    t0 = dt(2);
    dt = dt(1);
else,
    t0 = 0;
end

N = size(X,2);
if Reject50Hz,
    X = local_filter50(dt, X, N);
end

for icol=1:N,
    Xi = X(:,icol);
    Xi(Xi<Thr) = Thr;
    D = diff(Xi);
    iPeak = 1+find(D(1:end-1)>=0 & D(2:end)<0);
    iPeak = local_prune(iPeak, dt, Xi(iPeak), minISI);
    Spt{icol} = t0+dt*(iPeak-1);
    P {icol}= Xi(iPeak);
end


if isRow,
    X = X.';
end
%=========================
function  X = local_filter50(dt, X, N);
Time = Xaxis(X(:,1),dt/1e3); % time axis in s
Freq = 50; % Hz
Hum = exp(2*pi*i*Time*Freq);
for ii=1:N,
    Z = cosfit(dt, X(:,ii), Freq);
    Z = real(Z*Hum);
    X(:,ii) =  X(:,ii) - Z;
end


function ipeak = local_prune(ipeak, dt, Peak, minISI);
while 1,
    t_before = dt*[inf; diff(ipeak)]; % time between peak & previous one
    rp_before = [inf; diff(Peak)]; % relative peak height compared to previous one
    t_after = dt*[diff(ipeak); inf]; % ditto next one
    rp_after = [-rp_before(2:end); inf]; % relative peak height compared to next one
    % black list: peaks preceded by a recent, larger peak & peaks followed
    % by a soon, larger peak
    i_elim = find((t_before<minISI & rp_before<0) | (t_after<minISI & rp_after<0)); 
    if isempty(i_elim), break; end
    ipeak(i_elim) = [];
    Peak(i_elim) = [];
end % while 1



