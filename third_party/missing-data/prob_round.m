function y = prob_round(x)
% probabilistically round a numer.  ie, for x between 0 and 1,
% x is rounded up with prob. x.
f = floor(x);
y = f + ((x-f) > rand);
