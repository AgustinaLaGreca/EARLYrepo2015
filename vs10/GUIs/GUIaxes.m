function [h,AD]=GUIaxes(figh,Name);
% GUIaxes - handle of axes in GUI
%    GUIaxes(figh) returns the handle of the main AxesDisplay, which is 
%    the most recent AxesDisplay object named '@..' drawn in GUI with
%    handle figh.
%
%    GUIaxes(figh,'Foo') returns the handle to the AxisDisplay named Foo.
%
%    [h, AD] = GUIaxes(..) also returns the AxesDisplay object rendered as
%    axes object with handle h.
%
%    See also AxesDisplay.

if nargin<2,
    Name = ''; % main Axes display (see help text)
end

if isempty(Name),
    AD = getGUIdata(figh,'MainAxesDisplay');
else,
    AD = getGUIdata(figh,'AxesDisplay');
    AD = AD(Name);
end

h = AD.uiHandles.Axes;


