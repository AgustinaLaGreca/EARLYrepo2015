function reset(M)
% Messenger/reset - reset messenger object
%   reset(M) resets rendered Messenger M to its default values, The string
%   is set to M.TestLine (see Messenger); other property/value pairss are
%   taken from M.uicontrolProps. M may be an array of messengers.
%
%   See also Messenger, Messenger/report, Messenger/blank, GUImessage.

% delegate to report

for ii=1:length(M),
    report(M(ii),0);
end

