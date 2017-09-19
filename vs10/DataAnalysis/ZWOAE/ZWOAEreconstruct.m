function W = ZWOAEreconstruct(expID, recID, DPtype, Ref, plotflag)
%ZWOAEreconstruct -reconstruct ZWUIS stimulus based on ZWOAE recording
%
% W = ZWOAEreconstruct(expID, recID) takes the supall ZWOAEs (re. phase
% as send to speakers) and reconstructs the magnitude and phase of the
% zwuis components using wlinsolve. Weight factors are used
%
% W = ZWOAEreconstruct(expID, recID, DPtype) allows specification of a
% DP-type different than supall. Options are: 'near', 'far', 'supall',
% 'suplo' or 'suphi'.
%
% W = ZWOAEreconstruct(expID, recID, DPtype, Ref) allows specification of a
% phase reference. Default is 'stimulus', alternative is 'response' to obtain
% primary phases from the spectrum of the recorded signal.
%
% W = ZWOAEreconstruct(expID, recID, DPtype, Ref, plotflag) allows control
% over plotting. By default, no plot is generated. Two options for
% plotflag:
%   numeric: 1 (one) to plot, or 0 (zero) for no plot
%     <<OR>>
%   string: specifies plot arguments (e.g. 'r-p'); default is 'k-*'
%   This automatically assumes that plot has to be made.
%   Plots are added to current figure.
%
% Multiple recID are allowed for multiple reconstructions.
%
% Output W is a (array of) struct holding:
% expID:    experimental ID
% recID:    recording ID
% DPtype:   string specifying DPtype used in reconstruction
% Ref:      string specifying the phase reference used
% Fr:       ZWUIS frequencies [kHz]
% Mg:       'reconstructed' ZWUIS amplitudes [dB]
% Ph:       'reconstructed' ZWUIS phases [cycle]
%
%
% See also ZWOAEimport, ZWOAEmatrices, ZWOAEsubspec, wlinsolve
%

if nargin < 5 | isempty(plotflag), plotflag = 0; end
if nargin < 4 | isempty(Ref), Ref = 'stim'; end
if nargin < 3 | isempty(DPtype), DPtype = 'supall'; end

%...check of validity of Ref and DPtype is deferred to ZWOAEsubspec.m...

%--- use recursion for multiple recID ---
if numel(recID) > 1,
    for ii = 1:numel(recID),
        W(ii) = ZWOAEreconstruct(expID, recID(ii), DPtype, Ref, plotflag);
    end
    return;
end

%--- single recID from here ---

Q = ZWOAEsubspec(expID, recID, DPtype, Ref); %get the data
Ref = Q.Ref; DPtype = Q.subType; %to get 'full' names of these fields
M = ZWOAEmatrices(numel(Q.Fzwuis),['M' DPtype]); %get approp. ZWOAE matrix
M(:,end) = []; %remove last column; Fsingle is fixed

%unsort the DP-magnitudes, DP-phases and weigth factors to get them in order as
%obtained via ZWOAEmatrix*[Fzwuis; Fsingle]
Mg = Q.MG(Q.isortFreq);
Ph = Q.PH(Q.isortFreq);
Wght = Q.W(Q.isortFreq);

Mg = wlinsolve(Mg, M, Wght); %amplitude reconstruction for zwuis
Ph = wlinsolve(Ph, M, Wght); %phase reconstruction for zwuis
Mg = Mg(:); Ph = Ph(:); %force results into column arrays
Ph = Ph - round(mean(Ph)); %force mean-Ph to be on interval [-.5 5]

Fr = Q.Fzwuis; Fr = Fr(:); %get the ZWUIS frequencies [kHz]

%collect data in output struct W
if ~isreal(recID), recID = [Q.recID]; end
W = CollectInStruct(expID, recID, DPtype, Ref,'-', Fr, Mg, Ph);

%--- plotting ---
if ischar(plotflag) || isequal(plotflag, 1),
    if isnumeric(plotflag), plotflag = 'k-*'; end
    set(gcf,'units', 'normalized', 'position', [0.297 0.458 0.414 0.436]);
    subplot(2,1,1); box on;
    xplot(Fr, Mg-max(Mg),plotflag);
    ylabel('rel. amplitude [dB]');
    title(['expID: ' int2str(expID) '; recID: ' vector2str(recID),...
        '; PH-ref: ' Ref '; DPtype: ' DPtype]);

    subplot(2,1,2); box on; grid on;
    xplot(Fr, Ph, plotflag);
    xlabel('Fzwuis [kHz]'); ylabel('rel. phase [cycle]');
end

















