
%load the graph
load ./data/Dt_Sigmod.mat

%find the (AND) ceps between 'Jiawei_Han' and 'H._V._Jagadish',w/ budget 10
[Re] = CePS_Demo({'Jiawei_Han','H._V._Jagadish'},Wau,au_idx,10,2);