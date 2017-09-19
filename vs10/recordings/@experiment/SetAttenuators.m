function SetAttenuators(EXP, atten_dB, atten_dB_R, H);
% Experiment/SetAttenuators - set analog attenuators 
%     Experiment(EXP,[L R]) or SetAttenuators(EXP,L,R) sets the left- and
%     right-channel attenuators to L dB and R dB, respectively. The 
%     attenuator type is determined from the Experiment definition EXP. 
%     If no attenuators are present (EXP.Attenuators='-'), the 
%     only accepted attenuation values are 0 dB and 'max' (see below).
%     
%     SetAttenuators(EXP,X) sets both attenuators to X dB.
%
%     SetAttenuators(EXP, 'max') sets both attenuators to their maximum
%     value, which is determined from the CT.Hardware field.
%
%     SetAttenuators(EXP, L, 'max') sets the left attenuator to L dB, and 
%     the right one to its maximum value.
%     SetAttenuators(EXP, 'max', R) .. guess what.
%
%     See Experiment/edit, sys3PA5.


% [added, JV 13/8/2015]
if nargin<3,
    % only set both attenuators to X dB if both are present.
    switch lower(EXP.Audio.DAchannelsUsed)
        case 'both'
        % do nothing extra
        case 'left'
            if isnumeric(atten_dB), atten_dB = atten_dB(1); end;
            atten_dB = {atten_dB, nan};
        case 'right'
            if isnumeric(atten_dB), atten_dB = atten_dB(end); end;
            atten_dB = {nan, atten_dB};
    end
end

if nargin>2, 
    atten_dB = {atten_dB, atten_dB_R};
end
atten_dB = cellify(atten_dB);

if nargin>3, % overrule first arg
    clear EXP;
    EXP.Audio.Attenuators = H.AttenType;
end

% action is hardware specific
switch upper(EXP.Audio.Attenuators),
    case 'PA5', 
        
         atten_dB = local_max(atten_dB, 120); 
         
        % modified by EVE 14/08/2015 
        % right attenuator is not present in devicelist, exclude him by
        % putting NaN as value
        temp = sys3devicelist;
        attenuators = temp((strmatch('PA5', sys3devicelist))');
        if any(strcmp(attenuators,'PA5_2'))
             
        else
      atten_dB{2} = nan; 
        end

       
        sys3PA5(atten_dB{:});
        
        
    case '-', 
        atten_dB = local_max(atten_dB, 0);
        if ~isequal(0,unique([atten_dB{:}])),
            error('Cannot realize non-zero attenuation without attenuators.');
        end
    otherwise, error(['Unknown attenuator type ''' EXP.Audio.Attenuators '''.']);
end

function atten_dB = local_max(atten_dB, MaxAtten);
% replace 'max'
if isequal('max', atten_dB{1}),
    atten_dB{1} = MaxAtten;
end
if isequal('max', atten_dB{end}),
    atten_dB{end} = MaxAtten;
end










