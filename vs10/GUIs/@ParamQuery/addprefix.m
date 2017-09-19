function Q = width(Q, Prefix);
% ParamQuery/addprefix - add prefix to the names of paramquery objects.
%   Q=addprefix(Q,'Foo') adds prefix 'Foo' to the names of all the elements
%   of paramquery array Q.
%
%   See also ParamQuery, ParamPanel.

for ii=1:numel(Q),
    Q(ii).Name = [Prefix, Q(ii).Name];
end








