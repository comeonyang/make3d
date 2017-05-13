function pcut = addmcgormick(p)

pcut = p;

z = p.bilinears(:,1);
x = p.bilinears(:,2);
y = p.bilinears(:,3);

still_uncertain = find(abs(p.lb(z)-p.ub(z))>1e-8);
if ~isempty(still_uncertain)

    x_lb = p.lb(x);
    x_ub = p.ub(x);
    y_lb = p.lb(y);
    y_ub = p.ub(y);
    m = length(x);
    one = ones(m,1);
    general_vals =[x_lb.*y_lb one -y_lb -x_lb,x_ub.*y_ub one -y_ub -x_ub,-x_ub.*y_lb -one y_lb x_ub,-x_lb.*y_ub -one y_ub x_lb]';  %3
    general_cols = [one z+1 x+1 y+1 one z+1  x+1 y+1 one z+1 x+1 y+1 one z+1 x+1 y+1]';
    general_row = [1;1;1;1;2;2;2;2;3;3;3;3;4;4;4;4];
    quadratic_row = [1;1;1;2;2 ;2; 3; 3; 3];

    quadratic_cols =  [one  z+1 x+1 one z+1 x+1 one  z+1 x+1]';
    quadratic_vals = [-x_ub.*x_lb -one x_lb+x_ub x_lb.*y_lb one -y_lb-x_lb x_ub.*y_ub one -y_ub-x_ub]';

    m = 1+length(p.c);
    rows = [];
    cols = [];
    vals = [];
    nrow = 0;
%    dummy = ismembc(p.bilinears(:,1),p.r);
    for i =still_uncertain(:)'% 1:size(p.bilinears,1)
        x = p.bilinears(i,2);
        y = p.bilinears(i,3);
        if x~=y
            if 0%dummy(i)
                here = find(p.bilinears(i,1) == p.r);
                if p.s(here) > 0
                    rows = [rows;general_row(1:8,:)+nrow];
                    vals = [vals;general_vals(1:8,i)];
                    cols = [cols;general_cols(1:8,i)];
                    nrow = nrow + 2;
                elseif p.s(here)<0
                    rows = [rows;general_row(1:8,:)+nrow];
                    vals = [vals;general_vals(9:end,i)];
                    cols = [cols;general_cols(9:end,i)];
                    nrow = nrow + 2;
                end
            else
                rows = [rows;general_row+nrow];
                vals = [vals;general_vals(:,i)];
                cols = [cols;general_cols(:,i)];
                nrow = nrow + 4;
            end
        else
            if 0%dummy(i)
                %col = quadratic_cols(:,i);
                %val = quadratic_vals(:,i);
                here = find(p.bilinears(i,1) == p.r);
                if p.s(here) > 0
                    rows = [rows;quadratic_row(1:6,:)+nrow];
                    vals = [vals;quadratic_vals(4:end,i)];
                    cols = [cols;quadratic_cols(4:end,i)];
                    nrow = nrow + 2;
                elseif p.s(here)<0
                    rows = [rows;quadratic_row(1:3,:)+nrow];
                    vals = [vals;quadratic_vals(1:3,i)];
                    cols = [cols;quadratic_cols(1:3,i)];
                    nrow = nrow + 1;
                    %rows = [rows;general_row+nrow];
                    %vals = [vals;val];
                    %cols = [cols;col];
                    %nrow = nrow + 4;
                end
            else
                col = quadratic_cols(:,i);
                val = quadratic_vals(:,i);
                rows = [rows;quadratic_row+nrow];
                vals = [vals;val];
                cols = [cols;col];
                nrow = nrow + 3;
            end
        end
    end

    F_temp = sparse(rows,cols,vals,nrow,m);
    keep = find(~isinf(F_temp(:,1)));
    F_temp = F_temp(keep,:);
    pcut.F_struc = [F_temp;pcut.F_struc];
    pcut.K.l = pcut.K.l+size(F_temp,1);
end

