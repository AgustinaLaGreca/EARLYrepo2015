function y = SeqplayGrab_adc(G);
% SeqplayGrab_adc - set or get ADC-grabbing status of seqplay circuit
%   SeqplayGrab_adc returns the adc_grabstatus of the currently loaded 
%   seqplay circuit. When loading the circuit, the default status is 0, 
%   which means that recrding from ADC channels is disabled. 
%
%   SeqplayGrab_adc(1) or SeqplayGrab_adc('on') enables ADC recording.
%
%   SeqplayGrab_adc(0) or SeqplayGrab_adc('off') enables ADC recording.
%
%   See also SeqPlay, Seqplayinit, SeqplayInfo.


P = private_seqPlayInfo; % info shared by seqplayXXX

if nargin<1, % get
    y = sys3getpar('DoGrab_ADC', P.Dev);
else, % set
    if isequal('on', G), G = 1; 
    elseif  isequal('off', G), G = 0; 
    end
    sys3setpar(G, 'DoGrab_ADC', P.Dev);
end


