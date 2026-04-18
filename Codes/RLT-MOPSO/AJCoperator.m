function sop = AJCoperator(pop, rep, nPop, VarMin, VarMax, w, VarSize, nVar, it, MaxIt)
    empty_particle = pop(1);
    sop = repmat(empty_particle, 0, 1);
    
    D_space = mean(VarMax - VarMin); 
    R_max = D_space / (nPop^(1/nVar));
    R_min = 0.1 * R_max;
    Rt = R_max * (1 - it/MaxIt) + R_min * (it/MaxIt);
    
    all_positions = reshape([pop.Position], nVar, nPop)';
    distances = pdist2(all_positions, all_positions);
    
    rho = zeros(nPop, 1);
    for i = 1:nPop
        rho(i) = sum(distances(i, :) < Rt) - 1; 
    end
    
    rho_average = mean(rho);
    
    n = 1;
    for i = 1:nPop
        if rho(i) < rho_average
            outside_idx = find(distances(i, :) >= Rt);
            
            if ~isempty(outside_idx)
                [~, max_idx] = max(distances(i, outside_idx));
                j = outside_idx(max_idx);
                
                sop_particle = pop(i);
                
                r_i = rand(VarSize);
                sop_particle.Velocity = w * sop_particle.Velocity ...
                                      + r_i .* (pop(j).Position - sop_particle.Position);
                
                sop_particle.Position = sop_particle.Position + sop_particle.Velocity;
                
                for p = 1:nVar
                    sop_particle.Position(p) = max(sop_particle.Position(p), VarMin(p));
                    sop_particle.Position(p) = min(sop_particle.Position(p), VarMax(p));
                end
                
                sop_particle.Cost = MOP2(sop_particle.Position);
                sop(n, 1) = sop_particle;
                n = n + 1;
            end
        end
    end
end