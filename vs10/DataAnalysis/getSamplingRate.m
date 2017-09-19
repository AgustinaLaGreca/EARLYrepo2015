function [ fs ] = getSamplingRate( exp_name, ds_id )
%getSamplingRate Returns the sampling rate of a previously recorded dataset
%   Input arguments:
%       exp_name: a String containing the name of the experiment
%       ds_id: an integer corresponding to the id of the dataset(see first 
%               collumn of Databrowse)
%
%   Output arguments
%       fs: sampling frequency in Hertz
%
%   Created by: Adriaan Lambrechts 10/03/2017
%

ds = read(dataset,exp_name,ds_id);
fs = ds.Rec.RecordInstr.Fsam;

end

