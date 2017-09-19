function R=reportstimprogress(SP, figh, DT, Ax, Name, You, ds_name, varargin);
% reportstimprogress - construct reportstimprogress object
%    R=reportstimprogress(SP, figh, DT, Ax) returns a reportstimprogress 
%    object R whose task it is to report the progress of ongoing stimulus
%    presentation. R is a action object. Inputs:
%       SP: stimpresent object containing the stimulus presentation
%     figh: handle to the GUI to which the report is sent (see GUImessage)
%       DT: update interval in ms.
%       Ax: upload destination of R (see Dynamic).
%
%    reportstimprogress(SP, figh, DT, figh, 'Foo') is the same as
%    reportstimprogress(SP, figh, DT, GUIaccess(figh, 'Foo')).
%
%    R=reportstimprogress(SP, figh, DT, figh, 'Foo', 'You') makes the being
%    finished of R (see action/isready) dependent on the being finished of
%    the Action object named 'You' in the same figh GUI (see action/AfterYou).
%    Instead of a single You, a char string of names {'Me' 'You' 'Them'}
%    are also okay.
%
%    Type 'methods reportstimprogress -full' to see what can be done with R.
%
%  See also sortConditions, GUImessage, action, dynamic.

if nargin<5, Name = ''; end % 
if nargin<6, You = ''; end % 

Stutter = 0; % Stim Progress is dependant on wheiter the stutter function is on or not
if ~isempty(varargin)
   if  isfield(varargin{1},'Stutter')
       Stim = varargin{1};
       if strcmpi(Stim.Stutter,'on');
           Stutter = 1;
       end
   end
end

isvoid = 0;
status = 'initialized';
if nargin<1,
    [SP, figh, DT] = deal([]);
    isvoid = 1;
    status = 'void';
end
R = CollectInStruct(SP, figh,ds_name,Stutter);
R = class(R ,mfilename, action(status, DT));
% upload
if isempty(Name), ADR={Ax}; else, ADR = {Ax Name}; end
if ~isvoid, R=upload(R,ADR{:}); end
% if a You is specified, define finish dependency
if ~isempty(You), R=afteryou(R, Ax, You); end









