

function leader = SelectLeader(rep, beta)

    GI = [rep.GridIndex];
    
    OC = unique(GI);
    
    N = zeros(size(OC));
    
    for k = 1:numel(OC)
        N(k) = numel(find(GI == OC(k)));
    end
    
    P = exp(-beta*N);
    P = P/sum(P);
    
    sci = RouletteWheelSelection(P);
    
    sc = OC(sci);
    
    SCM = find(GI == sc);
    smi = randi([1 numel(SCM)]);
    
    sm = SCM(smi);
    
    leader = rep(sm);

end