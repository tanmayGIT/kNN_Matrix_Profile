<p align="center">
  For the details of the algorithm see, the this article </p>
<h1> Matrix Profile Based kNN Search over Large Time Series </h1> 
<h4> Tanmoy Mondal, Reza Akbarinia, Florent Masseglia <br/> </h4>
<h6> INRIA & LIRMM, Univ. Montpellier, France </h6>  <br/> <br/> <br/>



The interest of kNN search over 1NN search is shown in Fig. 5 and 6 of the article. To generate these plot, run the following code and follow the instructions mentioned inside these code :
```
-- Journal_Exp_1.m
-- Journal_Exp_2.m
```

The use case, mentioned in Fig. 8 can be obtained by running the code :
```
-- outliers_plot.m
```
Fig. 10 : The variation of computational time with the increasing number of kNN similarity search. The computational time analysis is shown for three different proposed approaches i.e. max based (refer to Section 4.2.2 ), sort based (refer to Section 4.2.1) and heap-max based (refer to Section 4.2.3 ) for kNN similarity search. The data used to create this Fig. is obtained by running the code :

```
-- Call_Self_Join_Exp_VariousAlgo_Tanmoy.m
```

To obtain the data related to computational time of our proposed algorithm, run the following code :
```
-- Call_Exp_kNN_Tanmoy.m
```

To obtain the data related to computational time of adapted STOMP (Yeh et al. ) algorithm, run the following code :
```
-- Call_Exp_kNN_Keogh_Adapted.m
```


The generate the plots, shown in Fig. 11, run the following code :
```
-- All_Results_5/read_plot_Graph.m
-- All_Results_5/read_plot_Graph_1.m
```

