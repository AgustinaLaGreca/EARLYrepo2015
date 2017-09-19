function Att = sys3PA5get(iatt);  
% sys3PA5get - get attenuator values of PA5s 
%  
%   Att = sys3PA5(x) returns the current settings of the attenuators in a
%   2-element row array Att. It is assumed that two PA5s are connected.
%   If no attenuators are connected, [] is returned.
%
%   sys3PA5(iatt) only returns the setting of attenuator(s) # iatt. Default
%   value of iatt is [1 2].
%  
%   See also sys3dev, sys3PA5.

% changed by EVE
    % if nargin<1, iatt=[1 2]; end

    % if isempty(strmatch('PA5', sys3devicelist)), % no attenuators connected
    %     Att = []; 
    %     return;
    % end

    % Att(1,1) = invoke(sys3dev('PA5_1'),'GetAtten');
    % Att(1,2) = invoke(sys3dev('PA5_2'),'GetAtten');
    % 
    % Att = Att(1,iatt);

% from here the new code
%check how many attenuattors are present in devicelist

if nargin<1, iatt=[1]; end

temp = sys3devicelist;
attenuators = temp((strmatch('PA5', sys3devicelist))');

if isempty(attenuators)
    Att = []; 
    return;
else
    for L = 1:length(attenuators)
        Att(1,L) = invoke(sys3dev(attenuators{L}),'GetAtten');
    end
end

Att = Att(1,iatt);
