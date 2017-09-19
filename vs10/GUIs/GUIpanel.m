function P = GUIpanel(Name, Title, varargin);
% GUIpanel - define, but do not draw, uipanel containing param queries.
%    P = ParamPanel(Name, Title, ...), where ... stands for properties of
%    the uipanel, creates a GUIpiece consisting of a uipanel. Use add to
%    fill the panel with paramqueries and/or other GUIpieces.
%    Default fontsize is 10; default fontweight is bold.
%
%    Examples
%       P = GUIpanel('CarrierFreq', 'Carrier frequency', 'fontsiz', 12);
%       P = GUIpanel('CarrierFreq', '', 'bordertype', 'none');
%    Note that properties may be abbreviated. The second example produces
%    an invisible uipanel, without a title and border.

GS = GUIsettings; GSP = GS.GUIpanel;
% select those fields of GSP and varargin that match uipanel properties
DefaultUIpanelprops = FullFieldnames(GSP,uipanelProperties, 'select'); 
ExplicitUIpanelProps = FullFieldnames(struct('Title', Title, varargin{:}),uipanelProperties);
% merge default & explicit uipanel properties. The latter take precedence
uipanelProps = structJoin(DefaultUIpanelprops, ExplicitUIpanelProps);

P = GUIpiece(Name, uipanelProps, [0.1 0.1], ...
    GS.GUIpanel.XYorigin ,GS.GUIpanel.XYlowrightMargin); 



