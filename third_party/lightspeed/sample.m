function x = sample(p, n)
%SAMPLE    Sample from categorical distribution.
% X = SAMPLE(P,N) returns a row vector of N integers, sampled according to the 
% probability distribution P (an array of numbers >= 0, whose sum is > 0).
% sum(P) does not have to be 1, but it must be > 0.
% X(i) ranges 1 to length(P).

% Written by Tom Minka
% (c) Microsoft Corporation. All rights reserved.

% Written by Min Sun July 25th
% adapt case when p contain Inf entries

if nargin < 2
  n = 1;
end

% process case for Inf entries in p ==============
DetX = find(p == Inf);
RandList = find(p ~= Inf);
RandList = RandList(:);
p = p( RandList);
% ================================================

if n < 10
  cdf = cumsum(p(:));
  if cdf(end) <= 0
    error('distribution is all zeros');
  end
  RandX = zeros(1,n);
  for i = 1:n
    RandX(i) = sum(cdf < rand*cdf(end)) + 1;
  end
else
  % large n method
  p = p(:);
  p = p/sum(p);
  h = sample_hist(p,n);
  RandX = zeros(1,n);
  i = [0 cumsum(h)'];
  % set RandX = [1 1 1 1 2 2 3 3 3 ... ]
  for k = 1:length(h)
    RandX((i(k)+1):i(k+1)) = k;
  end
  RandX = RandX(randperm(n));
end

% ================
x = [DetX(:); RandList(RandX)];
% ================
