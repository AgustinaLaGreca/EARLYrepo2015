function Lab=CondLabel(D, varargin);
% Dataset/CondLabel - char string describing ta stimulus condition
%    CondLabel(D, K) returns a label indicating the Kth stimulus
%    condition of datset D. The value and unit are taken from 
%    D.Stim.Presentation.X. Note that K counts the conditions "as defined", 
%    which is not necessarily the seem as the order "as played". 
%    If K is an array, CondLabel returns a cell array of strings.
%    K==0 means "all conditions."
%
%    CondLabel(d, K, 'Y') uses independent variable Y instead of the
%    default X, etc.
%
%   For arrays, a cell string is returned having the same size as D.
%
%   See also dataset/IDstring, Stimpresentx/CondLabel.

% delegate Stimpresentx/CondLabel
Lab = CondLabel(D.Stim.Presentation, varargin{:});




