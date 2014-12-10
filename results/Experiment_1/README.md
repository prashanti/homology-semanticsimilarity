## Experiment 1
>Goal -- Observe if the addition of WithHomology grouping classes improve the similarity scores of any orthologies gene pair comparisons.

>Grouping Classes added to the ontology : 
  * Without Homology grouping classes
    * has\_part some (PATO:0000001 and inheres\_in some (E))
  * With Homology grouping classes
    *  has\_part some (PATO:0000001 and inheres\_in some (E or homologous\_to E)) 
    *  has\_part some (PATO:0000001 and inheres\_in some (E))