## Experiment 2
>Goal -- Observe if the addition of WithHomology grouping classes improve the similarity scores of any orthologies gene pair comparisons.

>Grouping Classes added to the ontology : 
  * Without Homology grouping classes
    * has\_part some (PATO:0000001 and inheres\_in some (E))
  * With Homology grouping classes
    *  has\_part some (PATO:0000001 and inheres\_in some homologous\_to E) 

>Notation in SimilarityScores.tsv
 * Column 3 (Better Similarity With Homology	Similarity) 
   * 1 indicates similarity is better with homology grouping classes.
   * -1 indicates similarity is worse with homology grouping classes.
   * 0 indicates the scores are same in both cases.
    