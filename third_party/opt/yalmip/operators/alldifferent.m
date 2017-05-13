function varargout=alldifferent(varargin)

switch class(varargin{1})

    case 'double'
        x = varargin{1};
        x = sort(x(:));
        varargout{1} = all(diff(x) > 0);

    case 'sdpvar'
        varargout{1} = set(yalmip('addextendedvariable',mfilename,varargin{:}) == 1);

    case 'char'
        x = varargin{3};
        [nx,mx] = size(x);
        x = reshape(x,nx*mx,1);
       
        [M,m] = derivebounds(x);
        if 0 & all(M(1) == M) & all(m(1) == m) & (nx*mx == (M(1)-m(1) + 1))
            % No benefit seen in this model in this special case so we
            % don't use it to (reduce possibility of introducing a new bug) 
            d = binvar(nx*mx,nx*mx,'full');
            F = set(d*(m(1):M(1))' == x) + set(sum(d,1) == 1) + set(sum(d,2) == 1);
        else
            % Add constraint |x(i)-x(j)| > 1
            pairs = nchoosek(1:nx*mx,2);
            d = binvar(length(pairs),1);
            x1 = x(pairs(:,1));
            x2 = x(pairs(:,2));

            F = set(x1 - x2 > 1-(1+M(pairs(:,1))-m(pairs(:,2))).*d);
            F = F + set(x2 - x1 > 1-(1+M(pairs(:,2))-m(pairs(:,1))).*(1-d));
        end
        varargout{1} = F;
        varargout{2} = struct('convexity','milp','monotonicity','milp','definiteness','milp','extra','marker');
        varargout{3} = varargin{3};
end