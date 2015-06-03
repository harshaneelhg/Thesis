--This fold contains a matlab inplementation of Center-Piece Subgraph


--main-entrance: CePS-Demo

--You need GraphViz for show the resulting subgraph (http://www.graphviz.org/)

-- ./data/Dt_Sigmod.mat contains a sample data of authorship graph from SIGMOD conference

--Example
  [Re] = CePS_Demo({'Jiawei_Han','H._V._Jagadish'},Wau,au_idx,10,2)


--See the paper for the details of algorithm:
***********************************************************************
        Hanghang Tong, Christos Faloutsos
        Center-piece subgraphs: problem definition and fast solutions. 
        KDD 2006: 404-413
***********************************************************************


---LOG:

4/28/2008: 
   - fix the bug on k_softand computation (thanks a lot to Kensuke Onuma!)
   - add a makefile (try to type in 'make', you will see how to use this code)

