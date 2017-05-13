function [p,x_min,upper] = initializesolution(p);

x_min = zeros(length(p.c),1);
upper = inf;
if p.options.usex0
    x = p.x0;
    z = evaluate_nonlinear(p,x);
    residual = constraint_residuals(p,z);
    relaxed_feasible = all(residual(1:p.K.f)>=-p.options.bmibnb.eqtol) & all(residual(1+p.K.f:end)>=p.options.bmibnb.pdtol);
    if relaxed_feasible
        upper = p.f+p.c'*z+z'*p.Q*z;
        x_min = x;
    end
else
    p.x0 = zeros(length(p.c),1);
    x = p.x0;
    z = evaluate_nonlinear(p,x);
    residual = constraint_residuals(p,z);
    relaxed_feasible = all(residual(1:p.K.f)>=-p.options.bmibnb.eqtol) & all(residual(1+p.K.f:end)>=p.options.bmibnb.pdtol);
    if relaxed_feasible
        upper = p.f+p.c'*z+z'*p.Q*z;
        x_min = x;
    end   
end

