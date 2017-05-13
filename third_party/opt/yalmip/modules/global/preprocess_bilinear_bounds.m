function p = preprocess_bilinear_bounds(p)

% if ~isempty(p.bilinears)
%     quadratics = find(p.bilinears(:,2)  == p.bilinears(:,3));
%     for i = quadratics(:)'
%         if ismember(p.bilinears(i,2),p.binary_variables)
%             p.binary_variables = [p.binary_variables p.bilinears(i,1)];
%         elseif ismember(p.bilinears(i,2),p.integer_variables)
%             p.integer_variables = [p.integer_variables p.bilinears(i,1)];
%         end
%     end
% end

if ~isempty(p.integer_variables)
    for i = 1:size(p.bilinears,1)
        if ismember(p.bilinears(i,2),p.integer_variables)
            if ismember(p.bilinears(i,3),p.integer_variables)
                p.integer_variables = [p.integer_variables p.bilinears(i,1)];
            end
        end
    end
    if p.K.f > 0
        for i = 1:p.K.f
            if all(p.F_struc(i,:) == fix(p.F_struc(i,:)))
                involved = find(p.F_struc(i,2:end));
                if nnz(ismember(involved,p.integer_variables)) == length(involved)-1
                    p.integer_variables = [p.integer_variables involved];
                end
            end
        end
        p.integer_variables  = unique(p.integer_variables );
    end
end

if ~isempty(p.binary_variables)
    for i = 1:size(p.bilinears,1)
        if ismember(p.bilinears(i,2),p.binary_variables)
            if ismember(p.bilinears(i,3),p.binary_variables)
                p.binary_variables = [p.binary_variables p.bilinears(i,1)];
            end
        end
    end
    for i = 1:p.K.f
        if all(p.F_struc(i,:) == fix(p.F_struc(i,:)))
            involved = find(p.F_struc(i,2:end));
            if nnz(ismember(involved,p.binary_variables)) == length(involved)-1
                p.binary_variables = [p.binary_variables involved];
            end
        end
    end
    p.binary_variables  = unique(p.binary_variables );

end

if isempty(p.ub)
    p.ub = repmat(inf,length(p.c),1);
end
if isempty(p.lb)
    p.lb = repmat(-inf,length(p.c),1);
end
if ~isempty(p.F_struc)
    [lb,ub,used_rows] = findulb(p.F_struc,p.K);
    if ~isempty(used_rows)
        lower_defined = find(~isinf(lb));
        if ~isempty(lower_defined)
            p.lb(lower_defined) = max(p.lb(lower_defined),lb(lower_defined));
        end
        upper_defined = find(~isinf(ub));
        if ~isempty(upper_defined)
            p.ub(upper_defined) = min(p.ub(upper_defined),ub(upper_defined));
        end
        % Remove linear bounds
        used_rows = used_rows(find(~any(p.F_struc(p.K.f+used_rows,1+p.nonlinears),2)));
        not_used_rows = setdiff(1:p.K.l,used_rows);
        for i = 1:length(p.KCut.l)
            p.KCut.l(i) = find(not_used_rows==p.KCut.l(i) );
            p.originalModel.KCut.l(i) = find(not_used_rows==p.originalModel.KCut.l(i) );
        end
        if ~isempty(used_rows)
            p.F_struc(p.K.f+used_rows,:)=[];
            p.K.l = p.K.l - length(used_rows);
            used_rows = used_rows(used_rows <= p.originalModel.K.l);
            p.originalModel.F_struc(p.originalModel.K.f+used_rows,:)=[];
            p.originalModel.K.l = p.originalModel.K.l - length(used_rows);
        end
    end
end
p.lb(p.binary_variables) = max(0,p.lb(p.binary_variables));
p.ub(p.binary_variables) = min(1,p.ub(p.binary_variables));
p.lb(p.integer_variables) = ceil(p.lb(p.integer_variables));
p.ub(p.integer_variables) = floor(p.ub(p.integer_variables));

p = clean_bounds(p);
if ~isempty(p.bilinears)
    p = updatenonlinearbounds(p);
end
