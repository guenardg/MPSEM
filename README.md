# R Package MPSEM

Modelling Phylogenetic Signals Using Eigenvector Maps

Computation of Phylognetic Eigenvector Maps (**PEM**) and simulation of trait
evolution among phylogenetic trees or any directed graph representing
reticulated phylogenies.

**PEM** are sets of orthogonal basis vectors (eigenvectors) that are tailored to
represent trait evolving neutrally (i.e., with trait values changing smoothly
along the edges and nodes of the phylogenies) or non-neutrally (i.e., with trait
values shifting abruptly at the nodes of the phylogeny following adaptive shifts
related with changing environment or niche). Each basis (eigen) vector
represents a potential trait evolution pattern and linearly combining sets of
such patterns enables one to represent trait (or meta-trait) evolution
(phylogenetic modelling) or partial out phylogenetic variation when assessing
ecological hypothesis (e.g., testing functional trait correlation hypotheses).

**MPSEM** also features functions to simulate trait evolving neutrally or
non-neutrally, with specified optimal trait values (and shifts thereof) and
evolution rates, and along either "classical" phylogenetic trees (i.e., not
involving lateral gene transfer) or reticulated phylogenies represented by a
directed graph (i.e., involving lateral gene transfer through, e.g.,
hybridization). That functionality is useful to simulate the range of trait
values that can emerge from different trait evolution scenario.

Maintained by Guillaume Guénard -- Université de Montréal
