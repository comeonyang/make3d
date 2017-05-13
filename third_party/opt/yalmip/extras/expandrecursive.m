function [F_expand,failure,cause] = expandrecursive(variable,F_expand,extendedvariables,monomtable,where,level,options,method,extstruct,goal_vexity)
global DUDE_ITS_A_GP ALREADY_MODELLED REMOVE_THESE_IN_THE_END MARKER_VARIABLES OPERATOR_IN_POLYNOM
cause = '';
failure = 0;

if  ~alreadydone(getvariables(variable))
    % Request epigraph (or hypograph) for this variable
    [F_graph,properties,arguments,fcn] = model(variable,method,options,extstruct);

    % This is useful in MPT
    if ~isempty(F_graph)
       F_graph = tag(F_graph,['Expansion of ' fcn]);
    end
    
    % Bit of a hack (CPLEX sucks on some trivial problems, so we have to track
    % these simple models and fix them in the end)
    if ~isempty(properties)
        if isfield(properties,'extra')
            if strcmp(properties.extra,'marker')
                MARKER_VARIABLES = [MARKER_VARIABLES getvariables(variable)];
            end
        end
    end

    % Now check that this operator actually is convex/concave/milp
    if isequal(method,'graph')

        switch properties.convexity
            case {'convex','concave'}
                failure = ~isequal(properties.convexity,goal_vexity);
            case 'milp'
                if options.allowmilp
                    failure = 0;
                    method = 'milp';
                else
                    failure = 1;
                end
            case {'none','exact'}
                failure = 0;
        end

        if failure & options.allowmilp
            [F_graph,properties,arguments] = model(variable,'milp',options);            
            % This is useful in MPT
            if ~isempty(F_graph)
               F_graph = tag(F_graph,['Expansion of ' fcn]);
            end
            if isempty(F_graph)
                cause = [cause ', MILP model not available.'];
                return
            else
                failure = 0;
                method = 'milp';
            end
        elseif failure
            cause = ['Expected ' 'goal_vexity' ' function in '  where ' at level ' num2str(level)];
        end
    else
        if isempty(F_graph)
            cause = ['MILP model not available in ' where ' at level ' num2str(level)];
            failure = 1;
            return
        else
            failure = 0;
            method = 'milp';
        end
    end

    % We save variable models for future use
    if failure == 0
        done = save_model_expansion(method,F_graph,properties);
        if done
            return
        end
    end

    % Now we might have to recurse
    arguments = arguments(:);

    [ix,jx,kx] = find(monomtable(getvariables(variable),:));
    if ~isempty(jx) % Bug in 6.1
        if any(kx>1)
            OPERATOR_IN_POLYNOM = [OPERATOR_IN_POLYNOM extendedvariables(jx(find(kx>1)))];
        end
    end

    % Small pre-processing to speed-up large-scale problems (subsref sloooow)
    % with only linear arguments (such as norm(Ax-b) problems)
    if isa(arguments,'sdpvar')
        do_not_check_nonlinearity = is(arguments,'linear');
        if do_not_check_nonlinearity
            allvariables = getvariables(arguments);
            fullbasis = getbase(arguments);
            fullbasis = fullbasis(:,2:end);
            fullbasis_transpose = fullbasis';
        end
    else
        do_not_check_nonlinearity = 0;
    end

    j = 1;
    % Ok, here goes the actual recursive code  
    while j<=length(arguments) & ~failure

        if do_not_check_nonlinearity
%            usedvariables = find(fullbasis(j,:));
            usedvariables = find(fullbasis_transpose(:,j));
            expressionvariables = allvariables(usedvariables);
        else
            expression = arguments(j);
            expressionvariables = unique([depends(expression) getvariables(expression)]);
        end
        index_in_expression = find(ismembc(expressionvariables,extendedvariables));

        if ~isempty(index_in_expression)
            for i = index_in_expression
                if ~alreadydone(expressionvariables(i))
                    if do_not_check_nonlinearity
                      %  basis = fullbasis(j,expressionvariables(index_in_expression(i)));
                       % basis = fullbasis(j,usedvariables(index_in_expression(i)));
                        basis = fullbasis(j,usedvariables((i)));
                    else
                        basis = getbasematrix(expression,expressionvariables(i));
                    end

                    go_convex1 = (basis > 0) &  isequal(goal_vexity,'convex') & isequal(properties.monotonicity,'increasing');
                    go_convex2 = (basis <= 0) &  isequal(goal_vexity,'convex') & isequal(properties.monotonicity,'decreasing');
                    go_convex3 = (basis <= 0) &  isequal(goal_vexity,'concave') & isequal(properties.monotonicity,'increasing');
                    go_convex4 = (basis > 0) &  isequal(goal_vexity,'concave') & isequal(properties.monotonicity,'decreasing');

                    go_concave1 = (basis > 0) &  isequal(goal_vexity,'convex') & isequal(properties.monotonicity,'decreasing');
                    go_concave2 = (basis <= 0) &  isequal(goal_vexity,'convex') & isequal(properties.monotonicity,'increasing');
                    go_concave3 = (basis <= 0) &  isequal(goal_vexity,'concave') & isequal(properties.monotonicity,'decreasing');
                    go_concave4 = (basis > 0) &  isequal(goal_vexity,'concave') & isequal(properties.monotonicity,'increasing');

                    if go_convex1 | go_convex2 | go_convex3 | go_convex4
                        [F_expand,failure,cause] = expandrecursive(recover(expressionvariables(i)),F_expand,extendedvariables,monomtable,where,level+1,options,method,[],'convex');
                    elseif go_concave1 | go_concave2 | go_concave3 | go_concave4
                        [F_expand,failure,cause] = expandrecursive(recover(expressionvariables(i)),F_expand,extendedvariables,monomtable,where,level+1,options,method,[],'concave');
                    elseif isequal(properties.monotonicity,'exact')
                        [F_expand,failure,cause] = expandrecursive(recover(expressionvariables(i)),F_expand,extendedvariables,monomtable,where,level+1,options,method,[],goal_vexity);
                    else
                        if options.allowmilp
                            [F_expand,failure,cause] = expandrecursive(recover(expressionvariables(i)),F_expand,extendedvariables,monomtable,where,level+1,options,'milp',[],'milp');
                        else
                            failure = 1;
                            cause = ['Monotonicity required at ' where ' at level ' num2str(level)];
                        end
                    end

                end
            end
        end
        if ~do_not_check_nonlinearity  & ~DUDE_ITS_A_GP & ~options.expandbilinear & ~options.allownonconvex
            if isa(expression,'sdpvar')
                if degree(expression)~=1 &~is(expression,'sigmonial')
                    [Q,c,f,x,info] = quaddecomp(expression);
                    if info
                        failure = 1;
                        cause = ['Polynomial expression  at ' where ' at level ' num2str(level)];
                    else
                        eigv = real(eig(Q));
                        if ~all(diff(sign(eigv))==0)
                            failure = 1;
                            cause = ['Indefinite quadratic in ' where ' at level ' num2str(level)];
                        else
                            fail1 =  isequal(goal_vexity,'convex')  & all(eigv<=0) &  ~isequal(properties.monotonicity,'decreasing');
                            fail2 =  isequal(goal_vexity,'convex')  & all(eigv>0)  &  ~isequal(properties.monotonicity,'increasing');
                            fail3 =  isequal(goal_vexity,'concave') & all(eigv<=0) &  ~isequal(properties.monotonicity,'increasing');
                            fail4 =  isequal(goal_vexity,'concave') & all(eigv>0)  &  ~isequal(properties.monotonicity,'decreasing');

                            if fail1 | fail3
                                failure = 1;
                                cause = ['Concave quadratic encountered in ' where ' at level ' num2str(level)];
                            elseif fail2 | fail4
                                failure = 1;
                                cause = ['Convex quadratic encountered in ' where ' at level ' num2str(level)];
                            end
                        end
                    end
                end
            end
        end
        j = j+1;
    end
    F_expand = F_expand + F_graph;
end
