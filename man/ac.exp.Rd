% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/acdatabase_experiments.R
\name{ac.exp}
\alias{ac.exp}
\title{Experiment}
\description{
Experiment
}
\details{
An experiment entry (exp) is an environment containing some of:

\tabular{ll}{
name   \tab   \cr
id \tab  6 random capitals and numerics \cr
description \tab  \cr
comments \tab  \cr
results \tab See below \cr
}

The results field contains a list, each with the assay results from a single original titer table from a collaborating group. Each of these entries is itself a list, containing:
#'\tabular{ll}{
lab   \tab  Where the experiment was carried out \cr
assay \tab  HI or FRA \cr
conducted_by \tab Who carried out the assay \cr
date \tab \cr
comments \tab  \cr
file \tab Name of the original titer table file in databases repo  \cr
antigen_ids \tab  \cr
serum_ids \tab  \cr
titers \tab Titer table wihth '*' denoting a missing value \cr
meta \tab  \cr
}
Notably, the parent environment is always the empty environment; $.parent is the backbone-wise parent antigen entry; alteration$.parent is the origin antigen for a transplanted gene

It has classes c("acdatabase.ag", "acdatabase.entry", "environment")
}
