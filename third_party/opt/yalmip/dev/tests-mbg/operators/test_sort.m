function test_sort

x = sdpvar(4,1);
z = sdpvar(4,1);

[y,loc] = sort(x);

w = randn(4,1);

sol = solvesdp(set(-100 < x < 100)+set(z == y),norm(x-w,1));

mbg_asserttrue(sol.problem == 0);
mbg_asserttolequal(norm(sort(w)-double(z)),0,1e-4);


