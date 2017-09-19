function UCvphaseParam(keyword, varargin);
% UCvphasePARAM - parameter dialog of UCvphase

% use generic UCxxxParam script to handle generic stuff
callerfnc = mfilename; 
DialogName = 'VPHASEparam';
UCxxxParamSwitchboard; 


%========specific stuff is handles by local functions below===============

%================================================================
function localFillEdits(hh, Param);
% params are passed as struct; display them in appropriate edits
% without strict testing. This fcn should always be synchronized
% with UCvphase (it may be better to merge these two)
if nargin<2, Param = []; end;
% ---standard edits: Nbin Isub, timeWindow, Yunit, Ymax
[ds, hh, Param] = DataPlotFillEdits(hh, Param, 'VPHASE');
% ---wrap frequency for cycle histogram
% evaluate possible combinations of freq (car|mod|other) and ...
%  ... channel (beat|left|right|none) for current dataset
[ifc icc fstr cstr defaultstr] = ExploreFcycle(ds);
% initialize buttons accordingly
hfcyc = [hh.FcycButton, hh.ChanButton];
cycEdith = [hh.FcycEdit hh.FcycUnit]; % uicontrols whose visibilty depend on button 
UItwoButtons('init',hfcyc, ifc, icc, fstr, cstr, {cycEdith nan}, {[1 0 0 0] nan});
if isequal('auto', Param.FcycType),
   % set buttons to default settings obtained from exploreFcycle
   UItwoButtons('set', hfcyc, defaultstr{:});
else, % reproduce settings from last visit - do it in a foolproof way, though
   try,   
      UItwoButtons('set', hfcyc, Param.FcycType, Param.Chan);
      setstring(hh.FcycEdit, num2str(Param.Fcyc));
   catch, UItwoButtons('set', hfcyc, defaultstr{:});
   end
end
%---comp delay and Rayleigh crit
setstring(hh.CompDelayEdit, num2str(Param.compDelay));
setstring(hh.RayleighEdit, num2str(Param.RaleighCrit));
%---phase convention
set(hh.PhaseSignButton, 'userdata', Param.PhaseConv); MenuButtonMatch(hh.PhaseSignButton);
%---phase units
ii = strmatch(Param.PhaseUnit, {'cycle', 'rad', 'deg'});
set(hh.PhaseUnitButton, 'userdata', ii); MenuButtonMatch(hh.PhaseUnitButton);

function param = localReadEdits(hh);
if ~isstruct(hh), hh = getUIhandle(hh); end;
Nrep = getUIprop(hh.Root, 'Iam.ds.Nrep');
[dum, dum2, iRep] = IsubAndNbinCheck(hh, nan, Nrep);
UIinfo(''); % clear message window
param = [];
TimeWindow = TimeWindowCheck(hh);
if isempty(TimeWindow), return; end;
compDelay = UidoubleFromStr(hh.CompDelayEdit,1);
RaleighCrit = UidoubleFromStr(hh.RayleighEdit,1);
if ~checkNaNandInf([compDelay, RaleighCrit]), return; end;
% Fcycle
FcycType = getstring(hh.FcycButton);
Chan = getstring(hh.ChanButton);
visi = get(hh.FcycEdit, 'visible');
if isequal('on', visi), 
   Fcyc = UidoubleFromStr(hh.FcycEdit,1);
else, Fcyc = 100;
end
% phase convention
PhaseConv = UIintFromToggle(hh.PhaseSignButton); % 1|2 = lag|lead>0
%---phase units
PhaseUnit = get(hh.PhaseUnitButton, 'string'); % cycle|rad|deg
param = CollectInStruct(iRep, TimeWindow, RaleighCrit, compDelay, FcycType, Chan, ...
   Fcyc, PhaseConv, PhaseUnit);


