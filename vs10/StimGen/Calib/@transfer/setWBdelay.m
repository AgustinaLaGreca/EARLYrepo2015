function S = setWBdelay(S, tau);
% transfer/setWBdelay - set wideband delay of transfer function
%    S = setWBdelay(S, tau) sets the wideband delay of S to tau ms.
%    This is just a matter of bookkeeping; the phases of the complex
%    transfer fcn of T are changed accordingly.
%    If S is an array, the delayss of its elements are set to the
%    corresponding elements of tau.
%
%    S = setWBdelay(S, 'flat') sets the delay such that the remaining 
%    phase gradient in the transfer function vanishes. It does so by
%    fitting a straight line to the phase-freq curve. 
%
%    S = setWBdelay(S, 'wflat') uses a power-weighted fit insstead: the
%    weighting factors are propertional to the power of the transfer
%    function, Power=abs(H(freq)).^2.
%
%    S = setWBdelay(S, i*tau) sets the wideband delay of S to tau ms, but 
%    leaves the phases alone.
%
%
%    See Transfer, Transfer/measure, Transfer/getWBdelay.

if ischar(tau),
    if ~isequal('flat', lower(tau)) && ~isequal('wflat', lower(tau)),
        error('The only valid string values for tau are ''flat'' and ''wflat''.');
    end
    if isequal('flat', lower(tau)),
        W = 1+0*frequency(S); % constant weight of 1
    elseif isequal('wflat', lower(tau)),
        W = dB2P(magnitude(S)); % weights proportional to power
    end
    W = W/sum(W);
    for ii=1:numel(S),
        PP = wpolyfit(frequency(S(ii))/1e3, phase(S(ii)),W,1);
        Tau = -PP(1);
        S(ii) = setWBdelay(S(ii), getWBdelay(S(ii))+Tau);
    end
    return;
elseif isnumeric(tau) && isreal(i*tau), % just set stored WB_delay, no phase manipulation (see help)
    S.WB_delay_ms = imag(tau);
    return;
end

[S,tau] = SameSize(S,tau);

for ii=1:numel(S),
    dtau_s = 1e-3*(tau(ii)-S(ii).WB_delay_ms); % time shift in s
    Phasor = exp(2*pi*i*dtau_s*S(ii).Freq);
    S(ii).Ztrf = S(ii).Ztrf.*Phasor;
    S(ii).WB_delay_ms = tau(ii);
end
















