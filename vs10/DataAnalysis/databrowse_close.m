function databrowse_close( hObject,callbackdata, figh, Mode,gmess,gmessmode )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if ~isa(figh,'double')
    figh = figh.Number;
end;
GUImode(figh, Mode,gmess,gmessmode);
delete(hObject);
end

