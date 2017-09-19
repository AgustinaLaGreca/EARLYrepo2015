function [ Data ] = data( eventdata_obj )
%data gets the data from an eventdata object

    if ~isempty(eventdata_obj)
        Data = eventdata_obj.Data;
    end

end

