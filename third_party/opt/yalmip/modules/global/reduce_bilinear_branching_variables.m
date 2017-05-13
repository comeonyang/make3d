function p = reduce_bilinear_branching_variables(p);
if p.solver.lowersolver.objective.quadratic.convex
    % Setup quadratic
    Q_ = p.Q;
    for i = 1:size(p.bilinears,1)
        if p.c(p.bilinears(i,1))
            Q_(p.bilinears(i,2),p.bilinears(i,3)) = p.c(p.bilinears(i,1))/2;
            Q_(p.bilinears(i,2),p.bilinears(i,3)) = Q_(p.bilinears(i,3),p.bilinears(i,2))+p.c(p.bilinears(i,1))/2;
        end
    end
    if nnz(Q_)>0 & all(eig(full(Q_))>-1e-12)
        Used_in_F = find(any(p.F_struc(:,2:end),1));
        Used_in_F = intersect(Used_in_F,p.bilinears(:,1));
        p.branch_variables = [];
        for i = 1:size(p.bilinears,1)
            j = p.bilinears(i,1);
            if ismember(j,Used_in_F)
                p.branch_variables = [p.branch_variables p.bilinears(i,2:3)];
            end
        end
        p.branch_variables = unique( p.branch_variables);
    end
end