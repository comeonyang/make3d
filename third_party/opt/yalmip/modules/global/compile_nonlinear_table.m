function [linears,bilinears,nonlinears] = compile_nonlinear_table(p)
linears = find(sum(p.monomtable,2)==1);
nonlinears = find(~(sum(p.monomtable~=0,2)==1 & sum(p.monomtable,2)==1));
bilinears   = [];
for i = 1:length(nonlinears)
    z = find(p.monomtable(nonlinears(i),:));
    if length(z)==1
        bilinears = [bilinears;nonlinears(i) z z];
    else
        bilinears = [bilinears;nonlinears(i) z(1) z(2)];
    end
end

nonlinears = union(nonlinears,p.evalVariables);
linears = setdiff(linears,p.evalVariables);