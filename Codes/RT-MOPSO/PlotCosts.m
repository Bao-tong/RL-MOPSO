

function PlotCosts(pop, rep,Obj)


 rep_costs = [rep.Cost];    
if Obj==2
    plot(rep_costs(1, :), rep_costs(2, :), 'r*', 'MarkerSize', 8)   
    xlabel('1^{st} Objective');
    ylabel('2^{nd} Objective');
end
if Obj==3
    plot3(rep_costs(1, :), rep_costs(2, :), rep_costs(3, :),'r*', 'MarkerSize', 8)   
    xlabel('1^{st} Objective');
    ylabel('2^{nd} Objective');
    zlabel('3^{nd} Objective');
end
    grid on;
    
end