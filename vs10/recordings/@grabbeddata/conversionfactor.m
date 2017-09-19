function [CF, Unit] = conversionfactor(G);
% grabbeddata/conversionfactor - conversion factor grabbeddata
%   [F Unit]=conversionfactor(G) returns scale factor F and a string Unit
%   of the data in G. Methods of grabbeddata subclasses that return the
%   samples of G (e.g. adc_data/samples, adc_data/anmean, etc) return the
%   raw data, i.e., the Voltage V recorded by the AD conversion. To convert
%   to the physical units behind the recording - if archived -, multiply
%   the voltage by F. This will yield the physical units describe by Unit.
%   If no information is know about the physical units, conversionfactor
%   returns CF = 1 and Unit = 'V'.
%
%   See also grabbeddata/sourcedevice, grabbeddata/DataUnit.

[CF, V] = deal(1, 'V'); % defaults; see help text
S = sourcedevice(G);
if isfield(S, 'Settings'),
    if isfield(S.Settings, 'Unit'), Unit = S.Settings.Unit; end;
    if isfield(S.Settings, 'ScaleFactor'), 
        if ~isempty(strfind(S.ID.Type, 'Polytec_OFV5000')), % laser stuff: 
            CF = S.Settings.ScaleFactor;  % *range* is specified (e.g. 50 mm/s/V), not sensitivity
        else, % default: sensitivity is specified, e.g., 1 V/Pa
            CF = 1/S.Settings.ScaleFactor;
        end
    end;
end
if isnan(CF),
    CF = 1;
    Unit = 'V';
end









