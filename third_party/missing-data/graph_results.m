function o = graph_results(ifile, means, stds, ylabs,logscale)
% This will take two matrices with means and Standard deviations.
% We assume each row will be a single graph.
padpercent = .1;
% I'm going to pad the max and min values that I want visible in 
% the graph by this percent, so they don't touch the top and bottom.
if logscale
   if ~isempty(stds)
      stds_low = log(means - stds);
      stds_hi = log(means + stds);
   end
  means = log(means);
  ylabs = log(ylabs);
else
   if ~isempty(stds)
      stds_low = means - stds;
      stds_hi = means + stds;
   end
end

if isempty(stds)
   minval = min(min(ylabs),min(min(means)));
   maxval = max(max(means));
else
   minval = min(min(stds_low));
   maxval = max(max(stds_hi));
end
range = maxval - minval;

padded_minval = max(0,minval - padpercent*range);
padded_maxval = maxval + padpercent*range;

for i = 1:size(means,1)
   mean_row = means(i,:);
   if ~isempty(stds)
      std_low_row = stds_low(i,:);
      std_hi_row = stds_hi(i,:);
   else
      std_low_row = [];
      std_hi_row = [];
      end
  graph_idraw(ifile, mean_row, std_low_row, std_hi_row, padded_minval, padded_maxval);
end

label_graph(ifile, 1:size(means,2), 0, 1 + size(means,2), ylabs, padded_minval, padded_maxval);
o=1;
