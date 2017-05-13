function p = updatenonlinearbounds(p,changed_var,keepbest);
if ~isempty(p.bilinears)
    x = p.bilinears(:,2);
    y = p.bilinears(:,3);
    z = p.bilinears(:,1);
    x_lb = p.lb(x);
    x_ub = p.ub(x);
    y_lb = p.lb(y);
    y_ub = p.ub(y);
    bounds = [x_lb.*y_lb x_lb.*y_ub x_ub.*y_lb x_ub.*y_ub];
    new_lb = max([p.lb(z) min(bounds,[],2)],[],2);
    new_ub = min([p.ub(z) max(bounds,[],2)],[],2);
    % Avoid updating small bounds (numerical reasons)
    update = find(p.lb(z) < p.ub(z)-1e-4);
    p.lb(z(update)) = new_lb(update);
    p.ub(z(update)) = new_ub(update);
    
    p.lb(p.integer_variables) = fix(p.lb(p.integer_variables));
    p.ub(p.integer_variables) = fix(p.ub(p.integer_variables));
    p.lb(p.binary_variables) = fix(p.lb(p.binary_variables));
    p.ub(p.binary_variables) = fix(p.ub(p.binary_variables));

    quadratic_variables = p.bilinears(x==y,1);
    p.lb(quadratic_variables(p.lb(quadratic_variables)<0)) = 0;
end
