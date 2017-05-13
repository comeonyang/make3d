function x = evaluate_nonlinear(p,x)

if ~isempty(p.bilinears)
    x(p.bilinears(:,1)) = x(p.bilinears(:,2)).*x(p.bilinears(:,3));

    badones = find(sum(p.monomtable(p.bilinears(:,1),:),2)>2);
    if ~isempty(badones)
        redo = p.bilinears(badones,1);
        for i = 1:length(redo)
            these =  find(p.monomtable(redo(i),:));
            x(redo(i)) = prod(x(these).^((p.monomtable(redo(i),these))'));
        end
    end
else
    x = x(1:length(p.c));
    nonlinear = find(p.variabletype);
    x(nonlinear) = prod(repmat(x(:)',length(nonlinear),1).^p.monomtable(nonlinear,:),2);
end

for i = 1:length(p.evalMap)
    arguments = {p.evalMap{i}.fcn,x(p.evalMap{i}.variableIndex)};
    arguments = {arguments{:},p.evalMap{i}.arg{2:end-1}};
    x(p.evalVariables(i)) = feval(arguments{:});
end
