function test_alldifferent
n = 4;
x = sdpvar(n,1);
sol = solvesdp(set(1<x<n) + set(alldifferent(x)),sum(x))

mbg_asserttrue(sol.problem == 0)
mbg_asserttolequal(sort(double(x)),(1:n)', 1e-4);
