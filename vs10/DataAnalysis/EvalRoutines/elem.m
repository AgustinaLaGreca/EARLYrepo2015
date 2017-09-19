function element = elem(array, index)
%elem returns the element from array with index: 'index'
    element = [];
    if ~isempty(array)
        if index <= length(array)
           if isa(array,'cell')
               element = array{index}; 
           else
               element = array(index); 
           end
        end
    end
end

