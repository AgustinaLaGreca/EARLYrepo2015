function cleanup(D, varargout);
% grabbeddata/cleanup - cleanup temporary data associated with grabbeddata.
%   cleanup(D) removes any temporary data associated with G from disk.
%   For each concrete Grabbeddata subclass, Cleanup must be overloaded. It 
%   is an error to call the non-overloaded version.
%
%   See also adc_data/cleanup, Dataset/cleanup.

eval(IamAt);

error(['Cleanup is not yet overloaded for subclass ''' class(D) '''.' char(10) ...
    'It is an error to call the non-overloaded grabbeddata/cleanup.']);






