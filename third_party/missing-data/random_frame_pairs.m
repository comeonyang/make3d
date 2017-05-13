function randompairs = random_frame_pairs(n)
% Given n frames, pick all different pairs of frames, and randomize their order.
% This is returned in a 2xn matrix, where each column gives a pair of frames.
pairs = (all_ntuples(2,n))';
randompairs = pairs(1:2,randperm(size(pairs,2)));
