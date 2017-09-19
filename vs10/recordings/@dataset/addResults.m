function D = addResults(D,Freq_masker,dB_masker,thrs,thr_probe,cf,thr_data,SR)
% Dataset/addThr - add THR data to DataSet.
%
%   See Dataset, Experiment.
D.Status = 'complete';

D.Data(1).Freq_masker = Freq_masker;
D.Data(1).dB_masker = dB_masker;
D.Data(1).thrs = thrs;
D.Data(1).thr_probe = thr_probe;
D.Data(1).cf = cf;
D.Data(1).thr_data = thr_data;
D.Data(1).SR = SR;







