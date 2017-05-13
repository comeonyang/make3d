function F = SetupEvaluationVariable(varargin)

z = varargin{end};
X = varargin{3};
if isequal(getbase(X),[0 1])
    F = set([]);
else
    dX = double(X);
    if ~all(isnan(dX))
        assign(z,double(X));
    end
    F = set(X == z);
end