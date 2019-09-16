# FLIC
This repo hosts some of the code used for "Control of sugar and amino acid feeding via pharyngeal taste neurons" (Chen et al., currently in review).  

The first script (Simple Sorting Code.R) analyzed the raw electrical signal data of each well to generate feeding features described throughout this paper. A baseline electrical signal was calculated for each well by identifying the most commonly occurring electrical signal value throughout the assay period. Feeding events were defined as samples with electrical signal over 100 a.u. + baseline; tasting events were defined as samples with electrical signal over the baseline but less than 100 a.u. + baseline. Two consecutive events were distinguished if they occurred more than 1 second apart.
  Definition of features
-	The first feeding intensity: peak electrical signal of the first feeding event
-	The first feeding duration: time elapsed across the first feeding event
-	Interval between the first and second feeding events: time elapsed between the first and second feeding events
-	Time to the first feeding event: timestamp of the initiation of first feeding event


The second script (FLIC loading organization.R) annotated each trial with its corresponding genotype and tastant information by extracting information from accompanying text files using character matching/extracting functions. A sample text file has been included in the repo to demonstrate appropriate formatting. This is critical, because this script reads information via line numbers.

The third script (logistic regression.R) performed a multiple logistic regression analysis using a model trained on a set of 400 randomly selected and manually sorted usable and unusable observations from our larger dataset of over 3000 individual trials. The regression predicted a continuous variable between 0 and 1 as a function of the weighted features input into the model. The area under the Receiver operating characteristic (ROC) curve was close to 1.00 (0.94), indicating good predictive performance. The regression scored about 95% accuracy when testing against a test set of 400 randomly selected and manually sorted trials. The regression provided about 67% yield of usable trials when testing against our overall dataset.
