function [H,G]=threshold_curve(D, figh, P);
% dataset/threshold_curve - threshold curve, works only for THR datasets
%    See also dataviewparam, dataset/enableparamedit.

% handle the special case of parameter queries. Do this immediately to 
% avoid endless recursion with dataviewparam.
if isvoid(D) && isequal('params', figh),
    [H,G] = local_ParamGUI;
    return;
end

% open a new figure or use existing one?
if nargin<2 || isempty(figh),
    open_new = isempty(get(0,'CurrentFigure'));
    figh=gcf; 
    if ~isa(figh,'double')
        figh = double(figh);  % MH July 2019 conversion to double should not be needed now that isSingleHandle has been fixed
    end;
else,
    open_new = isSingleHandle(figh);
end

% parameters
if nargin<3, P = []; end
if isempty(P), % use default paremeter set for this dataviewer
    P = dataviewparam(mfilename); 
end

% delegate the real work to local fcn
H = local_threshold_curve(D, figh, open_new, P);

% enable parameter editing
enableparamedit(D, P, figh);



%============================================================
%============================================================
function data_struct = local_threshold_curve(D, figh, open_new, P);
P = struct(P);
P = P.Param;
% prepare plot
figure(figh); clf; ah = gca;

if strcmpi(P.FreqScale,'log')
    x = D.Stim.Freq;
    y = D.Data(1).Thr;
    semilogx(x,y);
else
    x = D.Stim.Freq;
    y = D.Data(1).Thr;
    plot(x,y);
end

freqs = x;
thrs=y;
data_struct.freqs = x;
data_struct.thrs = y;
save('C:\Early_StimDefLeuven\adaptive\tckl_thr\thr.mat','freqs','thrs')


xlabel('frequency (Hz)','fontsize',10);
ylabel('threshold (dB SPL)','fontsize',10);
title(IDstring(D, 'full'),'fontsize', 12, 'fontweight', 'bold', 'interpreter', 'none');
if open_new, placefig(figh, mfilename, D.Stim.GUIname); end

%=======================
function [T,G] = local_ParamGUI;
% return GUI for the params
%=====GUI==========
% the "Scale panel"
Scale = GUIpanel('Scale', 'Scale'); % , 'backgroundcolor', 'r'
% order of parameters for sorting the conditions. First indep var is "fastest moving" in display order
FreqScale = ParamQuery('FreqScale', 'Freq scale', '', {'Lin' 'Log'}, '', ...
    '', 1);
Scale = add(Scale, FreqScale);
%P_order = marginalize(P_order,[3 5]);
% GUI as a whole
G=GUIpiece([mfilename '_parameters'],[],[0 0],[10 4]);
G = add(G,Scale);
G=marginalize(G,[40 20]);
% list all parameters in a struct
T = VoidStruct('FreqScale');
















