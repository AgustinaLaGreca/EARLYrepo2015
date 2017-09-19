function Ti = inverse(T, D);
% transfer/inverse - inverse of transfer function
%    inverse(T) is the inverse transfer function of T.
%
%    inverse(T, 'Foo') also sets the Description of T.
%
%    See Transfer, Transfer/times.

Ti = transfer; % needed to make Ti a Transfer object

if nargin<2,
    D = ['inverse of ' T.Description];
end

Ti.Q_stim = T.Q_resp;
Ti.Q_resp = T.Q_stim;
Ti.Ref_stim = T.Ref_resp;
Ti.Ref_resp = T.Ref_stim;
Ti.dBref_stim = T.dBref_resp;
Ti.dBref_resp = T.dBref_stim;
Ti.Description = D;
Ti.CalibParam = T.CalibParam;
Ti.Freq = T.Freq;
Ti.Ztrf = 1./T.Ztrf;
Ti.WB_delay_ms = -T.WB_delay_ms;
    








