function [Comps, Weights, Approx] = mypca(W, Nmax);
% mypca - principal components analysis for waveforms
%    [Comps, Weights, Approx] = mypca(W, Nmax)
%     W is Nsam x Nwav waveform matrix; columns are individual waveforms.
%     Nmax is max order of analysis.
%     Outputs.
%         Comps: Nsam x Nmax matrix whose columns are the components
%       Weights: Nmax x Nwav weight matrix, such that Comps*Weights is
%                an approximation of W.
%        Approx: approximation of W, i.e., Comps*Weights with the column
%                means of W added.
[Nmax] = arginDefaults('Nmax', []);

[Nsam, Nwav] = size(W);
Nmax = replaceMatch(Nmax, [], Nwav)

colMean = mean(W,1); % means of coulmns which go lost in PCA
[Weights, Comps] = princomp(W);
Weights = Weights.';
% force normalization  for which squared sum of weights is unity for
% each reconstructed component. This is probably lready done by princomp,
% but its hepl text doesn't say say.
stdW = sqrt(sum(Weights.^2,2));
Weights = bsxfun(@times, Weights, 1./stdW);
Comps = bsxfun(@times, Comps, 1./stdW.');

Comps = Comps(:,1:Nmax);
Weights = Weights(1:Nmax,:);
Approx = Comps*Weights;
dsize(Approx, colMean);
Approx = Approx + repmat(colMean,[Nsam,1]);

















