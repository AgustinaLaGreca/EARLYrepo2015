function [Value, EditStr, Mess] = expr2str(hh,Q);
% expr2str = user defined function to interprete uidata
%   expr2str(hh,Q) extracts 'String' from handle hh, 
%   and returns it as a string Value, so that it can be treated as a
%   a regular Matlab expression later on as eval(Value)

Mess = '';
Value = get(hh,'string');
EditStr = Value;