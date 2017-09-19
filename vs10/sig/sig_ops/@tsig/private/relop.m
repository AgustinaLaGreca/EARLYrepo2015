function [P, Mess] = relop(S,T,op);
% tsig/relop - relational operators for tsig. Helper function
%   syntax:  [P, Mess] = relop(S,T,op)
%    See tsig/eq etc. Outisde=domain values are set to false by convention.

% work is delegated to private function binop
defS=0; defT=0; 
if all(haslogic(S)), defS=false; end
if all(haslogic(T)), defT=false; end
[P, Mess] = binop(S,T,op,defS,defT);
% binop returns keyword-style messages. Expand them.
switch Mess,
    case 'void', Mess = 'Void tsig objects may not be compared.';
    case 'numnonrowvector', Mess = 'Numerical values compared to tsig objects must be scalars or row vectors.';
    case 'nchan', Mess = 'Incompatible channel counts between sig objects or between tsig and numerical value.';
    case 'nontsig', Mess = 'Tsig objects can only be compared to numerical values or other tsig objects.';
    case 'fsam', Mess = 'Tsig object having different sample frequencies cannot be be compared.';
end





