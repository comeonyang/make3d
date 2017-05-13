function F = robustify_sdp_conic(F_xw,Zmodel,x,w)

if length(F_xw) == 0
    F = [];
    return;
else
    Ftemp = [];
    for i = 1:length(F_xw)
        Fi = sdpvar(F_xw(i));
        if degree(Fi(:),w) > 1
            [BilinearizeringConstraints,failure] = deriveBilinearizing(Fi,w);
            if failure
                error('Cannot get rid of nonlinear uncertainty in uncertain SDP')
            else
                Ftemp = Ftemp + BilinearizeringConstraints;
            end
        end
    end
    F_xw = F_xw + Ftemp;

    if any(Zmodel.K.q) | any(Zmodel.K.s)
        error('Only polytope uncertainty supported for uncertain SDPs');
    else
        % FIX : Assumes all uncertainty in all constraints
        K = Zmodel.K;
        A = -Zmodel.F_struc((1+K.f):(K.f + K.l),2:end);
        b =  Zmodel.F_struc((1+K.f):(K.f + K.l),1);

        try
            vertices = extreme(polytope(A,b))';
        catch
            disp('You probably need to install MPT (needed for vertex enumeration)')
            disp('http://control.ee.ethz.ch/~joloef/wiki/pmwiki.php?n=Solvers.MPT')
            error('MPT missing');
        end
        if K.f > 0
            Aeq = -Zmodel.F_struc(1:K.f,2:end);
            beq =  Zmodel.F_struc(1:K.f,1);
            feasible = sum(abs(Aeq*vertices - repmat(beq,1,size(vertices,2))),1) < 1e-6;
            vertices = vertices(:,feasible);
            if isempty(feasible)
                error('The uncertainty space is infeasible.')                
            end
        end

        F = set([]);
        for j = 1:length(F_xw)
            for i = 1:size(vertices,2)
                F = F + set(replace(F_xw(j),w,vertices(:,i)));
            end
        end
    end
end
