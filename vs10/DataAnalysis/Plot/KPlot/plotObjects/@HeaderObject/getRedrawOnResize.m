function doRedraw = getRedrawOnResize(HeaderObject)
doRedraw = HeaderObject.params.RedrawOnResize;
% 
% %% see if the handle exists
% hdlExists = 1;
% try
%     findobj(HeaderObject.hdl);
% catch
%     hdlExists = 0;
% end
% 
% %% return values
% if isequal(1, hdlExists)
%     doRedraw = HeaderObject.params.RedrawOnResize;
% else
%     doRedraw = 0;
% end
