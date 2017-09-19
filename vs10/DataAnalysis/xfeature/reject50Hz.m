function [X,Z] = reject50Hz(dt, X);
% reject50Hz - 50-Hz reject filter
%    [X,Z] = reject50Hz(dt, X);
%    Inputs:
%            dt: sample period in ms
%             X: recording trace. Matrices are treated as sets of rows.
%    Outputs:
%            X: filtered version of input X
%            Z: complex amplitude of 50-Hz subtracted sinusoid

[X, isRow] = TempColumnize(X);
N = size(X,2);
[X,Z] = local_filter50(dt, X, N);

if isRow,
    X = X.';
end
%=========================
function  [X,z] = local_filter50(dt, X, N);
Time = Xaxis(X(:,1),dt/1e3); % time axis in s
Freq = 50; % Hz
Hum = exp(2*pi*i*Time*Freq);
for ii=1:N,
    z(ii) = cosfit(dt, X(:,ii), Freq);
    Z = real(z(ii)*Hum);
    X(:,ii) =  X(:,ii) - Z;
end




