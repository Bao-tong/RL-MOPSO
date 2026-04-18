
clc;
clear;
close all;

global PROB_NAME; 

benchmarks = {'Viennet2', 'Viennet3'};

num_benchmarks = length(benchmarks);
TEST = 30;         
MaxIt = 100;       
Iter_Steps = floor(MaxIt / 10);


ResultsData = zeros(num_benchmarks * 2, TEST + 2);
IterData = zeros(num_benchmarks * 2, Iter_Steps);
rowNames_all = cell(num_benchmarks * 2, 1);

OpsData = zeros(num_benchmarks * 3, MaxIt);
rowNames_ops = cell(num_benchmarks * 3, 1);

for b = 1:num_benchmarks
    name = benchmarks{b};
    rowNames_all{2*b-1} = [name, '_HV'];
    rowNames_all{2*b}   = [name, '_IGD'];
    
    rowNames_ops{3*b-2} = [name, '_Gbest'];
    rowNames_ops{3*b-1} = [name, '_DPRT'];
    rowNames_ops{3*b}   = [name, '_IVT'];
end

colNames_res = cell(1, TEST + 2);
for i = 1:TEST
    colNames_res{i} = ['Test', num2str(i)];
end
colNames_res{TEST+1} = 'Mean';
colNames_res{TEST+2} = 'Std';

colNames_iter = cell(1, Iter_Steps);
for i = 1:Iter_Steps
    colNames_iter{i} = ['Iter', num2str(i * 10)];
end

colNames_ops = cell(1, MaxIt);
for i = 1:MaxIt
    colNames_ops{i} = ['Iter', num2str(i)];
end

for b = 1:num_benchmarks
    prob_name = benchmarks{b};
    PROB_NAME = prob_name;  
    
    fprintf('\n==================================================\n');
    fprintf(prob_name);
    fprintf('==================================================\n');
    
    if startsWith(prob_name, 'DTLZ')
        nVar = 22; 
        Obj = 3;
        VarMin = zeros(1, nVar);
        VarMax = ones(1, nVar);
    elseif startsWith(prob_name, 'WFG')
        nVar = 12; 
        Obj = 3;
        VarMin = zeros(1, nVar);
        VarMax = 2:2:2*nVar; 
    elseif startsWith(prob_name, 'Viennet')
        nVar = 2;  
        Obj = 3;
        VarMin = -4 * ones(1, nVar); 
        VarMax = 4 * ones(1, nVar);
    end
    VarSize = [1 nVar]; 
    
       PF_Path = ['true PF', prob_name, '.xlsx'];
    try
        [PF] = xlsread(PF_Path);
    catch
        error('can not read PF document', PF_Path);
    end
    
    CostFunction = @(x) MOP2(x);

    IIGD_current = zeros(1, TEST);
    HHV_current = zeros(1, TEST);
    iter30_IGD_current = zeros(1, Iter_Steps);
    iter30_HV_current = zeros(1, Iter_Steps);
    
    for test = 1:TEST
        fprintf(prob_name, test, TEST);
        
        nPop = 200; rows = 20; cols = 10; nRep = 200; w = 0.2; wdamp = 0.99;
        c1 = 1; c2 = 1.93; nGrid = 7; alpha = 0.1; beta = 2; gamma = 2;
        u_idx = 1; EPnum = 3; MPnum = 2;
        
        empty_particle.Position = []; empty_particle.Velocity = []; empty_particle.Cost = [];
        empty_particle.Best.Position = []; empty_particle.Best.Cost = []; empty_particle.IsDominated = [];
        empty_particle.Leader.Position = []; empty_particle.GridIndex = []; empty_particle.GridSubIndex = [];
        empty_particle.Normalization = []; empty_particle.Explore = []; empty_particle.Memory = [];
        empty_particle.repdis = []; empty_particle.topology = []; 
        pop = repmat(empty_particle, nPop, 1);

        for i = 1:nPop
            pop(i).Position = arrayfun(@(min, max) unifrnd(min, max), VarMin, VarMax);
            pop(i).Velocity = zeros(VarSize);
            pop(i).Cost = CostFunction(pop(i).Position);
            pop(i).Best.Position = pop(i).Position;
            pop(i).Best.Cost = pop(i).Cost;
        end

        pop = DetermineDomination(pop);
        rep = pop(~[pop.IsDominated]);
        Grid = CreateGrid(rep, nGrid, alpha);
        for i = 1:numel(rep)
            rep(i) = FindGridIndex(rep(i), Grid);
        end

        for it = 1:MaxIt
            pop = repdis(pop, nPop, rep, Obj);
            
            if it > 1
                sop = AJCoperator(pop, rep, nPop, VarMin, VarMax, w, VarSize, nVar, it, MaxIt);
                pop = [pop; sop; rep];
            end
            
            if numel(pop) > nPop
                Extra = numel(pop) - nPop;
                for e = 1:Extra
                    pop = DeleteOneRepMemebr(pop, gamma);
                end
            end
            
            pop = Resortexplore(pop, nPop, Obj, EPnum);
            pop = Resortmemory(pop, nPop, Obj, MPnum);
            
            if test == TEST
                count_op = zeros(1, 3);
            end
            
            for i = 1:nPop
                rt = rand;
                if rt >= 0 && rt < 0.333
                    pop(i).Leader.Position = Leaderparticle(pop, i, nPop, rows, cols);       
                    pop(i).Velocity = w * pop(i).Velocity ...
                        + c1 * rand(VarSize) .* (pop(i).Best.Position - pop(i).Position) ...
                        + c2 * rand(VarSize) .* (pop(i).Leader.Position - pop(i).Position);
                    pop(i).topology = 1;
                    
                elseif rt >= 0.333 && rt < 0.666
                    pop(i).Best.Position = Memoryparticle(pop, i, nPop);
                    pop(i).Leader.Position = Leaderparticle1(pop, i, nPop);       
                    pop(i).Velocity = w * pop(i).Velocity ...
                        + c1 * rand(VarSize) .* (pop(i).Best.Position - pop(i).Position) ...
                        + c2 * rand(VarSize) .* (pop(i).Leader.Position - pop(i).Position);
                    pop(i).topology = 2;
                    
                else
                    leader = SelectLeader(rep, beta);       
                    pop(i).Velocity = w * pop(i).Velocity ...
                        + c1 * rand(VarSize) .* (pop(i).Best.Position - pop(i).Position) ...
                        + c2 * rand(VarSize) .* (leader.Position - pop(i).Position);      
                    pop(i).topology = 3;
                end
                
                if test == TEST
                    count_op(pop(i).topology) = count_op(pop(i).topology) + 1;
                end
                
                pop(i).Position = pop(i).Position + pop(i).Velocity;
                for p = 1:nVar
                    pop(i).Position(p) = max(pop(i).Position(p), VarMin(p));
                    pop(i).Position(p) = min(pop(i).Position(p), VarMax(p));
                end
                pop(i).Cost = CostFunction(pop(i).Position);
                
                if Dominates(pop(i), pop(i).Best)
                    pop(i).Best.Position = pop(i).Position;
                    pop(i).Best.Cost = pop(i).Cost;
                elseif ~Dominates(pop(i).Best, pop(i))
                    if rand < 0.5
                        pop(i).Best.Position = pop(i).Position;
                        pop(i).Best.Cost = pop(i).Cost;
                    end
                end
            end
            
            if test == TEST
                row_idx_ops = 3 * b - 2;
                OpsData(row_idx_ops, it)     = count_op(3) / nPop; % Gbest
                OpsData(row_idx_ops + 1, it) = count_op(2) / nPop; % DPRT
                OpsData(row_idx_ops + 2, it) = count_op(1) / nPop; % IVT
            end
            
            rep = [rep; pop(~[pop.IsDominated])]; 
            rep = DetermineDomination(rep);  
            rep = rep(~[rep.IsDominated]);   
            Grid = CreateGrid(rep, nGrid, alpha);
            for i = 1:numel(rep)
                rep(i) = FindGridIndex(rep(i), Grid);
            end    
            if numel(rep) > nRep        
                Extra = numel(rep) - nRep;
                for e = 1:Extra
                    rep = DeleteOneRepMemebr(rep, gamma);
                end       
            end
            w = w * wdamp;
            
            if mod(it, 10) == 0
                PopObj = reshape([rep.Cost], Obj, [])';
                if test == TEST 
                    iter30_IGD_current(u_idx) = IGD(PopObj, PF);
                    iter30_HV_current(u_idx) = HV(PopObj, PF);
                end
                u_idx = u_idx + 1;
            end
        end 

        PopObj_final = reshape([rep.Cost], Obj, [])';
        IIGD_current(test) = IGD(PopObj_final, PF);
        HHV_current(test) = HV(PopObj_final, PF);
        
        row_idx_hv = 2 * b - 1;
        row_idx_igd = 2 * b;
        
        ResultsData(row_idx_hv, test) = HHV_current(test);
        ResultsData(row_idx_hv, TEST+1) = mean(HHV_current(1:test));
        ResultsData(row_idx_hv, TEST+2) = std(HHV_current(1:test));
        
        ResultsData(row_idx_igd, test) = IIGD_current(test);
        ResultsData(row_idx_igd, TEST+1) = mean(IIGD_current(1:test));
        ResultsData(row_idx_igd, TEST+2) = std(IIGD_current(1:test));
        
        if test == TEST
            IterData(row_idx_hv, :) = iter30_HV_current;
            IterData(row_idx_igd, :) = iter30_IGD_current;
        end
        
        writetable(array2table(ResultsData, 'VariableNames', colNames_res, 'RowNames', rowNames_all), ...
                   'RT-MOPSO_Results.csv', 'WriteRowNames', true);
                   
        writetable(array2table(IterData, 'VariableNames', colNames_iter, 'RowNames', rowNames_all), ...
                   'RT-MOPSO_Iterations.csv', 'WriteRowNames', true);
                   
        writetable(array2table(OpsData, 'VariableNames', colNames_ops, 'RowNames', rowNames_ops), ...
                   'RT-MOPSO_Operators.csv', 'WriteRowNames', true);
                   
    end 
end 
