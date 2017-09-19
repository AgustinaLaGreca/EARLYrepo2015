function cleanup(D, varargout);
% adc_data/cleanup - cleanup temporary data associated with adc_data object.
%   cleanup(D) removes any temporary data associated with D from disk.
%   These are the disk dumps of the hoard objects used for realtime data
%   storage (see hoard/cat). Cleanup is called when the final data are not
%   saved.
%
%   See also grabbeddata/cleanup, dataset/cleanup.

eval(IamAt);

sam = D.Data.Samples;
if sam.NchunkSaved>0,
    delete([sam.Filename '*.bin']);
end











