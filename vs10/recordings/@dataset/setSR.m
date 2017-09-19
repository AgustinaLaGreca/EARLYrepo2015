function [ ds ] = setSR( ds, SR )
%setSR Sets the SR struct as a field of ds.Data
if ~isempty(ds)
       ds.Data(1).SR = SR;
end


end

