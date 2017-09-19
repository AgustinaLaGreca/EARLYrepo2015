function PPS = ZWOAEplotStruct(Prop);
% ZWOAEplotStruct - plot styles for ZWOAE data
%
%   ZWOAEplotStruct returns a struct with fields:
%
%   Spec, Floor, Prim, DPnear, DPfar, DPsuplo, DPsuphi, DPsupall, DPall
%
%   Each of these fields is a struct containing plot properties conventionally
%   used for the coressponding elements of ZWOAE spectra.
%   These fields may be passed to plot commands.
%
%   ZWOAEplotStruct(Prop) returns only the struct for the specified Prop.
%   Prop must be a string giving one of the fieldnames listed above.
%
% Examples:
%
%   1)
%   PS = ZWOAEplotStruct;
%   plot(Freq_near, Phase_near, PS.DPnear);
%
%   2)
%   line(AllFreqs, AllAmps, ZWOAEplotStruct('Spec'));
%
%
%   See also ZWOAEmatrices, KeywordMatch
%

%either zero or one input variable for function
error(nargchk(0, 1, nargin, 'struct'));

%list of all plot properties specified within function.
allProps = {'Spec', 'Floor', 'Prim', 'DPnear', 'DPfar', 'DPsupall',...
    'DPsuplo', 'DPsuphi', 'DPall' 'Zwuis' 'Zwuis2' 'Zwuis3'};

PPS = [];
% ------define LINE properties-------
% -----------------------------------

% stimulus-related spectral components
PS.Spec = struct('color',[1 1 1]*.7,'LineWidth',2);
% noise floor
PS.Floor = struct('color',[1 1 1]*.2,'LineWidth',2);
% Primaries
PS.Prim = struct('color','r', 'markerfacecolor','r', 'marker', 'o',...
    'linestyle', 'none', 'markersize', 4, 'linewidth',0.5);
% DPs near zwuis group
PS.DPnear = struct('color','m', 'markerfacecolor','m', 'marker', 'p',...
    'linestyle', 'none', 'linewidth',0.5);
% DPs near lonely primary
PS.DPfar = struct('color',[0 0 0.7], 'markerfacecolor',[0 0 0.7],...
    'marker', 'p', 'linestyle', 'none', 'linewidth',0.5);
% suppression
PS.DPsupall = struct('color','g', 'markerfacecolor','g', 'marker', 's',...
    'linestyle', 'none', 'markersize', 4, 'linewidth',0.5);
% LOWsuppression
PS.DPsuplo = struct('color',[0 0.5 0.2], 'markerfacecolor',[0 0.5 0.2],...
    'marker', 'v', 'linestyle', 'none', 'markersize', 5, 'linewidth',0.5);
% HIGHsuppression
PS.DPsuphi = struct('color',[0.2 0.5 0], 'markerfacecolor',[0.2 0.5 0],...
    'marker', '^', 'linestyle', 'none', 'markersize', 5, 'linewidth',0.5);
% All types of DPs
PS.DPall = struct('color',[0 0 0], 'markerfacecolor',[0 0 0],...
    'marker', 's', 'linestyle', 'none', 'markersize', 3, 'linewidth',0.5);
% zwuis itself (just like prim)
PS.Zwuis = struct('color','r', 'markerfacecolor','r', 'marker', 'o',...
    'linestyle', 'none', 'markersize', 4, 'linewidth',0.5);
% zwuis itself (just like prim)
PS.Zwuis = struct('color','r', 'markerfacecolor','r', 'marker', 'o',...
    'linestyle', 'none', 'markersize', 4, 'linewidth',0.5);
PS.Zwuis2 = struct('color','r', 'markerfacecolor','r', 'marker', 's',...
    'linestyle', 'none', 'markersize', 4, 'linewidth',0.5);
PS.Zwuis3 = struct('color','r', 'markerfacecolor','r', 'marker', 'v',...
    'linestyle', 'none', 'markersize', 4, 'linewidth',0.5);


% -----------------------------------
% -----------------------------------

%Return single plot properties in PPS when requested.
if nargin>0,
    switch lower(Prop),
        case 'dpzwuis', Prop = 'zwuis';
        case 'dpzwuis2', Prop = 'zwuis2';
        case 'dpzwuis3', Prop = 'zwuis3';
    end
    [Prop, Mess] = keywordMatch(Prop, allProps,'ZWOAE plot property:');
    error(Mess);
    PPS = PS.(Prop);
else
    PPS = PS; %default, return all plot properties in PPS
end


