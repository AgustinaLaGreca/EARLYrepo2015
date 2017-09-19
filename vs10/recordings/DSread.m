function DS = DSread(varargin);
%   DSread - read Early dataset from file.
%      DSread('RG09211', 6) reads dataset #6 from experiment RG09211.
%      
%      DSread('RG09211', inf) or DSread('RG09211') reads the most recently
%      saved dataset from experiment RG09211.
%      
%      DSread(E, iRec) reads recording # iRec from experiment E. Here E is
%      an Experiment object, e.g., E = current(experiment) works.
%      
%      The real work is delegated to dataset/read.
%
%      See also Dataset, Experiment.

DS = read(dataset, varargin{:});




