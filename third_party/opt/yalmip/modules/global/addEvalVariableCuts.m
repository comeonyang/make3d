function pcut = addEvalVariableCuts(p)

temp = p.F_struc(1:p.K.f,:);
p.F_struc(1:p.K.f,:) = [];
pcut = p;
for i = 1:length(p.evalMap)
    y = p.evalVariables(i);
    x = p.evalMap{i}.variableIndex;
    xL = p.lb(x);
    xU = p.ub(x);

    switch p.evalMap{i}.fcn
        case 'exp'
            [Ax, Ay, b] = bound_exp(xL,xU);
        otherwise
            Ax = zeros(0,1);
            Ay = zeros(0,1);
            b = zeros(0,1);
    end
    F_structemp = zeros(size(b,1),length(p.c)+1);
    F_structemp(:,1+y) = -Ay;
    F_structemp(:,1+x) = -Ax;
    F_structemp(:,1) = b;
    pcut.F_struc = [F_structemp; pcut.F_struc];
    pcut.K.l = pcut.K.l + length(b);

end
pcut.F_struc = [temp;pcut.F_struc];
