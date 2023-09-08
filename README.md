# rodent-tail-morphology
data and code for "The evolution of rodent tail morphology"


This is the data, tree file, and code necessary to run the 11 main models (and a version of the tail length model with an interaction term) presented in the manuscript "The evolution of rodent tail morphology" by Sheard, Skinner, and Caro.

To run suborder analyses, simply take the relevant subset after loading in the data; do double-check that all fixed effects are necessary, and restructure the prior specification accordingly.

Note that phangorn's method of forcing trees to be ultrametric currently (9/2023: R 4.3.1) un-roots them; users will need to either use the "force.ultrametric" command from phytools ( tree <- force.ultrametric(tree, method=c("extend")) ) or use an older version of R (version 4.1.3 works fine).

The contents of the .csv file are identical to those of the .xlsx file published with the manuscript. We also include the metadata (found as an additional tab in the .xlsx file published with the manuscript) here as a separate .csv file.

The trees are taken from Upham et al. 2019 (https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000494); they have been trimmed to include only rodents, though the code will trim them further.
