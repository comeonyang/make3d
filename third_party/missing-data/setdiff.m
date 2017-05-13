function c = setdiff(a,b)
c = [];
for i = a
  if ~ member(i,b)
     c = [c,i];
  end
end