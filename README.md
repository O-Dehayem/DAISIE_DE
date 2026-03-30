# DAISIE-DE: Backward Likelihood Reformulation of DAISIE

DAISIE-DE is a reformulation of the DAISIE (Dynamic Assembly of Island Biotas through Speciation, Immigration and Extinction) model for island diversification.
This repository contains the R code and C++ code for estimating colonization and diversification parameters using a
**backward (tip-to-root) likelihood approach**.

## Overview

The original DAISIE model estimates colonization, speciation, and extinction rates on islands using phylogenetic trees, assuming all species diversify at the same rates.
Its likelihood is calculated forward in time from the crown age to the present, conditional on extant species. While effective, this forward-time formulation is not readily
extensible to trait-dependent diversification models.

DAISIE-DE addresses this by reformulating the likelihood in **reverse time**: probabilities are assigned to species at the tips and updated backward from the present 
to the root. This backward formulation is analogous to approaches used in SSE (State-Speciation and Extinction) models and lays the groundwork for future 
trait-dependent extensions.

### When to Use


DAISIE-DE is designed to estimate **colonization, cladogenesis, anagenesis, and extinction rates** from phylogenetic trees of entire island communities. 
