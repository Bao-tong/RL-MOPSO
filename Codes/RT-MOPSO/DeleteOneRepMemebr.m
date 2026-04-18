function rep = DeleteOneRepMemebr(rep, gamma)

    GI = [rep.GridIndex];
    
    OC = unique(GI);
    
    N = zeros(size(OC));
    
    for k = 1:numel(OC)
        N(k) = numel(find(GI == OC(k)));
    end
    
 
    max_val = max(gamma * N);
    P = exp(gamma * N - max_val);
    P = P / sum(P);
    
    sci = RouletteWheelSelection(P);
    
    if isempty(sci)
        sci = randi(numel(OC));
    end
    
    sc = OC(sci);
    
    SCM = find(GI == sc);
    
    smi = randi([1 numel(SCM)]);
    
    sm = SCM(smi);
    
    rep(sm) = [];

end