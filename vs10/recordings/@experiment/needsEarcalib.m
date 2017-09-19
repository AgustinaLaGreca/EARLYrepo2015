function N=needsEarcalib(E);
% experiment/needsEarcalib - true if ear calibration needs to be done
%   needsEarcalib(E) returns True if the type of calibraion of E is "In
%   situ" and there are no ear calibration data yet.
%
%   See also Experiment/Edit, Experiment/calib_index.


N = isequal('In situ', E.Audio.CalibrationType) && isequal(0, calib_index(E));


