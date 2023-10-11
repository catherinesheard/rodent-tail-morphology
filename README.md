# rodent-tail-morphology

Main manuscript: Sheard, C., Skinner, N., & Caro, T. The evolution of rodent tail morphology. The American Naturalist.

Contact: Catherine Sheard, catherine.sheard@abdn.ac.uk

%%%%

This repository contains the data, tree file, and code necessary to run the 11 main models (and a version of the tail length model with an interaction term) presented in the manuscript, as well as two worked-through examples of how to do suborder analyses.

Within the "main models" folder, the "Autotomy.R" script contains the most annotations; users are recommended to begin there.

This code has been verified to run on R version 4.1.3, along with packages MCMCglmm version 2.35 and phangorn 2.10.0. If you would like to run this code on a more recent version of R (e.g., version 4.3.1), you will also need the package phytools (verified to work with version 1.9.16) and use the command "tree <- force.ultrametric(tree, method=c("extend"))" instead of the command "tree<-nnls.tree(cophenetic(tree),tree,rooted=TRUE)".

The contents of the .csv file are identical to those of the .xlsx file published with the manuscript. We also include the metadata (found as an additional tab in the .xlsx file published with the manuscript) below.

The trees are taken from Upham et al. 2019 (https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000494); they have been trimmed to include only rodents, though the code will trim them further.

%%%%

Metadata for the RodentData.csv:

Column name -> Explanation

CES-ID -> Internal ID

MSW05_Binomial -> Binomial name

TreeTipNamePartial -> Upham et al. tree tip name to match to databases that index based on genus-species

UphamTreeName.full -> Upham et al. tree tip name within the file itself

Suborder -> Rodent suborder

White -> Whether the tail is white (1) or not (0)

White-tips -> Whether the tail has a white tip (1) or not (0)

Black-tips -> Whether the tail has a black tip (1) or not (0)

Contrasting -> Whether the tail has a contrasting tip (1) or not (0)

Tufts -> Whether the tail has a fluffy tuft (1) or not (0)

Naked -> Whether the tail is naked (1) or not (0)

Fluffy -> Whether the tail is fluffy (1) or not (0)

Tail_length -> Tail length, as a percentage of  body length

Prehen-Bin -> Whether the tail is prehensile (1) or not (0)

Fatty -> Whether the tail has fatty deposits (1) or not (0)

Autotomy -> Whether the tail is autotomous (1) or not (0)

Mean_Annual_Temp -> Average range mean annual temperature

Noc -> Whether the species is nocturnal (1) or not (0)

Di -> Whether the species is diurnal (1) or not (0)

Crep -> Whether the species is crepuscular (1) or not (0)

SocialGrpSize -> Typical social group size variable, with 1 = solitary species, 2 = groups of two to five individuals, 3 = six to 25 individuals, 4 = 26-100 individuals, and 5 = >100 individuals

Litter_size -> Average litter size

Hib -> Whether the species hibernates (1) or not (0)

A -> Whether the species is arboreal (1) or not (0)

F -> Whether the species is fossorial (1) or not (0)

Q -> Whether the species is aquatic (1) or not (0)

T -> Whether the species is terrestrial (1) or not (0)

shade_score -> Habitat shade score, as defined in the manuscript

Openness-Score -> Habitat openness score, as defined in the manuscript

Mass -> Body mass, from PHYLACINE

PrecipVar -> Average range precipitation variation, as described in the manuscript

SnakeThreat	-> The 'snake predation threat' variable, as described in the manuscript

OwlThreat -> The 'owl predation threat' variable, as described in the manuscript

RaptorThreatNonOwl -> The non-owl 'raptor predation threat' variable, as described in the manuscript




