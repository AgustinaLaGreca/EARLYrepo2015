function SetAttenuators(CT, atten_dB, atten_dB_R, Hardware);
% stimcontext/SetAttenuators - set analog attenuators 
%     SetAttenuators(CT,[L R]) or SetAttenuators(CT,L,R) sets the left- and
%     right-channel attenuators to L dB and R dB, respectively. The 
%     attenuator type is determined from the Hardware field of stimulus 
%     context CT. If no attenuators are present (attenType='NONE'), the 
%     only accepted attenuation values are 0 dB and 'max' (see below).
%     
%     SetAttenuators(CT,X) sets both attenuators to X dB.
%
%     SetAttenuators(CT, 'max') sets both attenuators to their maximum
%     value, which is determined from the CT.Hardware field.
%
%     SetAttenuators(CT, L, 'max') sets the left attenuator to L dB, and 
%     the right one to its maximum value.
%     SetAttenuators(CT, 'max', R) .. guess what.
%
%     SetAttenuators(CT, L, R, H) uses hardware settings from struct H and
%     ignores CT.
%
%     See stimcontext, sys3PA5.

if nargin>2, 
    atten_dB = {atten_dB, atten_dB_R};
end
if nargin>3,
    clear CT;
    CT.Hardware = Hardware;
end
atten_dB = cellify(atten_dB);

MaxAtten = CT.Hardware.MaxAtten;
% replace 'max'
if isequal('max', atten_dB{1}),
    atten_dB{1} = MaxAtten;
end
if isequal('max', atten_dB{end}),
    atten_dB{end} = MaxAtten;
end
% call hardware-specific function
switch upper(CT.Hardware.AttenType),
    case 'PA5', sys3PA5(atten_dB{:});
    case 'NONE', 
        if ~isequal(0,unique([atten_dB{:}])),
            error('Cannot realize non-zero attenuation without attenuators.');
        end
    otherwise, error(['Unknown attenuator type ''' CT.Hardware.AttenType '''.']);
end











