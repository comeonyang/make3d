function [g,geq,dg,dgeq] = fmincon_congp(x,dummy1,dummy2,dummy4)

persistent prob
persistent A

if nargin == 4
    prob = dummy4;
    A = spalloc(max(prob.map),size(prob.A,1),0);
    for i = 1:max(prob.map)
        ind = find(prob.map==i);
        A(i,ind) = prob.b(ind)';
    end
    return
end

g = A*exp(prob.A*x);

ind = find(g<1e-2);
if ~isempty(ind)
    g(ind) = exp(log(1e-2)+(g(ind)-1e-2)/1e-2);
end
g = log(g);

if length(prob.h) > 0
    geq = log(prob.h) + prob.G*x;
    dgeq = prob.G;
else
    geq = [];
    dgeq = [];
end
% Should be correct, but it fails for some problems (test_gp_5)
% dg = [];
% z = prob.A*x;
% for i = 1:1:length(g)    
% %     zz = prob.A*x + log(A(i,:))';
%     % dg = [dg ((1/(sum(exp(zz))))*exp(zz)'*prob.A)'];
%     b = A(i,:)';
%     expz = exp(z);
%     expz(isinf(expz)) = 1e5;
%     dg = [dg ((1/(sum(b.*expz)))*(b.*expz)'*prob.A)'];    
% end
% %dg = dg';
% full(dg)