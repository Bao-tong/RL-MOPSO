# Reinforcement learning-guided multi-objective particle swarm optimization with adaptive topology
1. Parameter_tuning: parameter tuning experiments for twelve multi-objective algorithms based on SMAC3. 
   Norm_Neg_HV = - ( HV(Obtained_Pareto_Front) / HV(True_Pareto_Front) ); 
   Norm_Neg_HV stands for Normalized Negative Hypervolume. Norm_Neg_HV is the primary metric used in this repository to evaluate and compare the performance of multi-objective algorithms based on training sets.
2. Codes: full codes of six algorithms based on topological structure: MOPSO-Gbest, MOPSO-DPRT, MOPSO-Von-Neumann, RT-MOPSO and RL-MOPSO. Full codes of MOPSOCD, SMPSO, dMOPSO, MOEADD, RVEA and PESAII can be found at    PlatEMO (https://github.com/BIMK/PlatEMO)
3. Comparison_results:
   (1) IGD and HV tets results of twelve algorithms on 17 benchmark instances (30 times); 
   (1) Wilcoxon_results for twelve algorithms; 
   (2) Friedman_results for twelve algorithms; 
   (3) Holm_results for twelve algorithms. 
4. Case_study_results: six algorithms based on topological structures results based on two cases.
5. Effectiveness of the Q-learning: three operators(Gbest. DPRT and IVT)' selection probabilities of RL-MOPSO and RT-MOPSO based on DTLZ2, WFG8, and Viennet3
