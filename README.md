# Accounting for reporting delays in real-time phylodynamic analyses with preferential sampling

This repository contains the code to replicate the analyses conducted in the paper *Accounting for reporting delays in real-time phylodynamic analyses with preferential sampling*.

## Dependencies

For 00-obtain-genalogies.qmd:

- [R](https://www.r-project.org/)
- [mafft](https://mafft.cbrc.jp/alignment/software/)
- [BEAST](https://beast.community/)

For 01-full-analysis-from-paper.qmd:

- [R](https://www.r-project.org/)
- [INLA](https://www.r-inla.org/download-install)
- [phylodyn2](https://github.com/CatalinaMedina/phylodyn2) can install in R with `devtools::install_github("CatalinaMedina/phylodyn2")`

## Navigation

```
├── analysis                              
|   ├── data 
|       ├── Washington-genealogies        <- Genealogies for the full, observed,
|                                            and truncated data scenarios
|       ├── delays-for-simulations        <- Reporting delays from real data for
|                                            use in simulations
|       ├── gisaid_supplemental_table.pdf <- Provided by GISAID for sequences
|                                            used in Washington analysis
|       ├── us-states-covid-ny-times.csv  <- covid case data used to compare
|                                            results from Washington analysis
|       ├── washington-historic-delays.csv<- sampling and reporting dates for 
|                                            most recent month of samples
|       └── washington-sample-details.csv <- Sample summaries for each data
|                                            scenario
│   ├── 00-obtain-genealogies.qmd         <- Code to produce genealogies
│   ├── 01-full-analysis-from-paper.qmd   <- Code for analysis in manuscript 
│   ├── my-plotting-functions.R           <- My plot helper functions used in
|                                            01-full-analysis.qmd
│   └── simulation-functions.R            <- My sim helper functions used in 
|                                            01-full-analysis.qmd
└──     
```