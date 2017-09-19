function I=irec(D);
% Dataset/irec - recording number within experiment of dataset
%   irec(D) or D.irec returns the recording number of dataset D, i.e. its
%   index within the experiment. An array of indices is returned when D is
%   a dataset array.
%
%   By convention, void datasets have irec=0.
%
%   See also Dataset, dataset/irec_pooled.

for ii=1:numel(D),
    d = D(ii);
    if ispooled(d), 
        d = members(d);
        d = d(1);
    end
    if isvoid(d),
        I(ii) = 0;
    else,
        I(ii) = d.ID.iDataset;
    end
end
I = reshape(I, size(D));












