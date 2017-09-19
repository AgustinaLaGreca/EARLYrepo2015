function [W] = setAnnotation(W,index,Label, BeginOfBaseline)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    if nargin == 3
        W.Annotations.Label{index} = Label;
    elseif nargin == 4

        W.Annotations.Label{index} = Label;
        W.Annotations.BeginOfBaseline = BeginOfBaseline;
    else 
        dbstack;
        error('wrong number of input arguments');
    end


end

