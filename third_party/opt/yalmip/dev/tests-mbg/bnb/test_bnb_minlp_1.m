function minlp

randn('seed',12345);
rand('seed',12345);

x = sdpvar(5,1);
A = randn(15,5);
b = rand(15,1)*10;

obj = sum(x) + sum((x-3).^4);
sol = solvesdp(set(A*x < b) + set(integer(x)),obj,sdpsettings('bnb.solver','fmincon','warning',0))

mbg_asserttolequal(sol.problem,0);
