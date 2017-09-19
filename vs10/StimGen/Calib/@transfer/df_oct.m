function df = df_oct(T);
% transfer/df_oct - frequency spacing of complex transfer function 
%    df_oct(T) returns the frequency spacing of T in octaves. 
%
%    Note: It is assumed that the frequencies have logarithmic spacing, 
%    which is true when they have been created using transfer/measure and
%    transfer/patch.
%
%    See Transfer, Transfer/measure, Transfer/Frequency, Transfer/Phase.

df = log2(max(T.Freq)/min(T.Freq))/numel(T.Freq);





