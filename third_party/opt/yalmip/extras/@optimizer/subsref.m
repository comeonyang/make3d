function [u,problem,message] = subsref(self,subs)

if isequal(subs.type,'()')
    
    if length(subs.subs{1} == self.n)
        
    elseif size(subs.subs{1},2) == self.n
        subs.subs{1} = subs.subs{1}';
    else
        error('Input argument has wrong size');
    end        
    u = [];
    for i = 1:size(subs.subs{1},2)
        self.model.F_struc(1:self.n,1) = subs.subs{1}(:,i);
        eval(['output = ' self.model.solver.call '(self.model);']);
        u = [u output.Primal(self.map)];
        problem = output.problem;
        message = yalmiperror(output.problem);
    end
else
    error('Reference type not supported')
end