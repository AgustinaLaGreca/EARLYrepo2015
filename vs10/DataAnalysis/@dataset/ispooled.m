function ip = ispooled(DS);
% Dataset/ispooled - true for pooled dataset
%    ispooled(DS) returns true of DS is a pooled dataset.
%    Arrays result in an logical array [ispooled(DS(1)) ispooled(DS(2)) .. ]
%
%    See also dataset.

ip = [];
for ii=1:numel(DS),
    ip(ii) = isa(DS(ii), 'pooled_dataset');
end
ip = reshape(ip, size(DS));







