function [model,changed] = bilinearize(model)

% Assume we don't do anything
changed = 0;

% Are there really any non-quadratic terms?
if any(model.variabletype > 2)
    % Bugger...
    changed = 1;

    % Find a higher order term
    first_polynomial = find(model.variabletype == 3);
    model = fixbounds(model);
    first_polynomial = first_polynomial(1);
    powers = model.monomtable(first_polynomial,:);
    if nnz(powers) == 1
        model = univariate_bilinearize(model,first_polynomial,powers);
    else
        model = multivariable_bilinearize(model,first_polynomial,powers);
    end
%     
%     % Find inverses etc
%     first_sigmonial = find(model.variabletype == 4);
%     if ~isempty(first_sigmonial)
%         first_sigmonial = first_sigmonial(1);
%         powers = model.monomtable(first_sigmonial,:);
%         if any(powers ~= fix(powers)
%             error('model class not supported')
%         else
%             powers_new(powers>0) = 0;
%             powers_new(powers>0) = 0;
%         end
%     end
end

function   model = univariate_bilinearize(model,first_polynomial,powers);

% Fix initial?
fix_initials = ~isempty(model.x0);

% variable^power
variable = find(powers);
p1 = floor(powers(variable)/2);
p2 = ceil(powers(variable)/2);

powers_1 = powers;powers_1(variable) = p1;
powers_2 = powers;powers_2(variable) = p2;

% Only recursive if power>4
switch p1+p2
    case 3
        [model,index2] = findoradd(model,powers_2);
        % Now define new variable y, replace x^3 with x*y, and add
        % constraint y == x^2
        model.monomtable(end+1,end+1) = 1;
        model.variabletype(end+1) = 0;
        model.monomtable(first_polynomial,variable) = 1;
        model.monomtable(first_polynomial,end) = 1;
        model.variabletype(first_polynomial) = 1;
        model.F_struc = [zeros(1,size(model.F_struc,2));model.F_struc];
        model.K.f = model.K.f + 1;
        model.F_struc(1,end+1) = 1;
        model.F_struc(1,1+index2) = -1;
        model.c(end+1) = 0;
        model.Q(end+1,end+1) = 0;
        if fix_initials
            model.x0(end+1) = initial(model.x0,powers_2);
%            model.x0(end+1) = 0;
        end
        bound = powerbound(model.lb,model.ub,powers_2);
        model.lb(end+1) = -bound;
        model.ub(end+1) = bound;

    case 4
        [model,index2] = findoradd(model,powers_2);
        model.monomtable(end+1,end+1) = 1;
        model.variabletype(end+1) = 0;
        model.monomtable(first_polynomial,variable) = 0;
        model.monomtable(first_polynomial,end) = 2;
        model.variabletype(first_polynomial) = 2;
        model.F_struc = [zeros(1,size(model.F_struc,2));model.F_struc];
        model.K.f = model.K.f + 1;
        model.F_struc(1,end+1) = 1;
        model.F_struc(1,1+index2) = -1;
        model.c(end+1) = 0;
        if fix_initials
            model.x0(end+1) = initial(model.x0,powers_2);
           % model.x0(end+1) = 0;
        end
        model.Q(end+1,end+1) = 0;
        bound = powerbound(model.lb,model.ub,powers_2);
        model.lb(end+1) = 0;
        model.ub(end+1) = bound;
    otherwise
        [model,index1] = findoradd(model,powers_1);
        [model,index2] = findoradd(model,powers_2);
        model.monomtable(end+1,end+1) = 1;
        model.monomtable(end+1,end+1) = 1;
        model.variabletype(end+1) = 0;
        model.variabletype(end+1) = 0;
        model.monomtable(first_polynomial,variable) = 0;
        model.monomtable(first_polynomial,end) = 1;
        model.monomtable(first_polynomial,end-1) = 1;
        model.variabletype(first_polynomial) = 1;

        model.F_struc = [zeros(1,size(model.F_struc,2));model.F_struc];
        model.K.f = model.K.f + 1;
        model.F_struc(1,end+1) = 1;
        model.F_struc(1,1+index1) = -1;
        model.F_struc = [zeros(1,size(model.F_struc,2));model.F_struc];
        model.K.f = model.K.f + 1;
        model.F_struc(1,end+1) = 1;
        model.F_struc(1,1+index2) = -1;

        model.c(end+1) = 0;
        model.Q(end+1,end+1) = 0;
        model.c(end+1) = 0;
        model.Q(end+1,end+1) = 0;

        if fix_initials
            x0 = model.x0;
            model.x0(end+1) = initial(x0,powers_1);
            model.x0(end+1) = initial(x0,powers_2);
%            model.x0(end+1) = 0;
%            model.x0(end+1) = 0;
        end
        bound1 = powerbound(model.lb,model.ub,powers_1);
        bound2 = powerbound(model.lb,model.ub,powers_2),
        
        model.lb(end+1) = -bound1;
        model.ub(end+1) = bound1;
        model.lb(end+1) = -bound2;
        model.ub(end+1) = bound2;
end
model = bilinearize(model);




function model = multivariable_bilinearize(model,first_polynomial,powers);

% Fix initial?
fix_initials = ~isempty(model.x0);

variables = find(powers);
mid = floor(length(variables)/2);
variables_1 = variables(1:mid);
variables_2 = variables(mid+1:end);
powers_1 = powers;
powers_2 = powers;
powers_1(variables_2) = 0;
powers_2(variables_1) = 0;

[model,index1] = findoradd(model,powers_1);
if sum(powers_1)>1
    model.monomtable(end+1,end+1) = 1;
    pos1 = size(model.monomtable,2);
    model.variabletype(end+1) = 0;
    bound = powerbound(model.lb,model.ub,powers_1);
    model.lb(end+1) = -bound;
    model.ub(end+1) = bound;
    if fix_initials;model.x0(end+1) = initial(model.x0,powers_1);end
end

[model,index2] = findoradd(model,powers_2);
if sum(powers_2)>1
    model.monomtable(end+1,end+1) = 1;
    pos2 = size(model.monomtable,2);
    model.variabletype(end+1) = 0;
    bound = powerbound(model.lb,model.ub,powers_2);
    model.lb(end+1) = -bound;
    model.ub(end+1) = bound;    
    if fix_initials;model.x0(end+1) = initial(model.x0,powers_2);end
end

model.monomtable(first_polynomial,:) = 0;
if sum(powers_1)>1
    model.monomtable(first_polynomial,pos1) = 1;
else
    model.monomtable(first_polynomial,variables_1) = 1;
end
if sum(powers_2)>1
    model.monomtable(first_polynomial,pos2) = 1;
else
    model.monomtable(first_polynomial,variables_2) = 1;
end
model.variabletype(first_polynomial) = 1;

if sum(powers_1)>1
    model.F_struc = [zeros(1,size(model.F_struc,2));model.F_struc];
    model.K.f = model.K.f + 1;
    model.F_struc(1,end+1) = 1;
    model.F_struc(1,1+index1) = -1;
    model.c(end+1) = 0;
    model.Q(end+1,end+1) = 0;
end

if sum(powers_2)>1
    model.F_struc = [zeros(1,size(model.F_struc,2));model.F_struc];
    model.K.f = model.K.f + 1;
    model.F_struc(1,end+1) = 1;
    model.F_struc(1,1+index2) = -1;
    model.c(end+1) = 0;
    model.Q(end+1,end+1) = 0;
end
model = bilinearize(model);


function [model,index1] = findoradd(model,powers_1);

if length(powers_1) < size(model.monomtable,2)
    powers_1(size(model.monomtable,1)) = 0;
end

index1 = findrows(model.monomtable,powers_1);

if isempty(index1)
    model.monomtable = [model.monomtable;powers_1];
    model.monomtable(end,end+1) = 0;
    index1 = size(model.monomtable,1);
    model.c(end+1) = 0;
    model.Q(end+1,end+1) = 0;
    model.F_struc(1,end+1) = 0;
    bound = powerbound(model.lb,model.ub,powers_1);
    model.lb(end+1) = -bound;
    model.ub(end+1) = bound;    
    if ~isempty(model.x0)
        model.x0(end+1) = 0;
    end
    switch sum(powers_1)
        case 1
            model.variabletype(end+1) = 0;
        case 2
            model.variabletype(end+1) = 2;
        otherwise
            model.variabletype(end+1) = 3;
    end
end



function model = fixbounds(model);
polynomials = find(model.variabletype > 0);
LU = max([abs(model.lb) abs(model.ub)],[],2);
for i = 1:length(polynomials)
    j = polynomials(i);
    if j<=length(model.lb)
    vars  = find(model.monomtable(j,:));
    bound = 1;
    for k = 1:length(vars)
        bound = bound * LU(vars(k))^model.monomtable(j,vars(k));
    end
    model.lb(j) = max(model.lb(j),-bound);
    model.ub(j) = min(model.ub(j),bound);
    end
end

function bound = powerbound(lb,ub,powers)
LU = max([abs(lb) abs(ub)],[],2);
vars  = find(powers);
bound = 1;
for k = 1:length(vars)
    bound = bound * LU(vars(k))^powers(vars(k));
end

function z = initial(x0,powers)
z = 1;
vars  = find(powers);
for k = 1:length(vars)
    z = z * x0(vars(k))^powers(vars(k));
end




 