#' PK_example dataset
#' 
#' Plasma drug concentration over time for 1 individual.
#' @name PK.example
#' @docType data
#'
#'@usage data(PK.example)
#'
#'@details
#'  The data consists of a dataframe named PK.example where the plasma drug concentrations (mg/L, CONC column) over time (hours, time column)  for one individual (ID column) are stored.
#'  The individual receives a dose of 150mg (amt column) at time 0. 
#'  
#'@examples
#'\dontrun{
#'  data(PK.example)
#'  print(PK.example)
#'}
#'@keywords datasets
NULL

#' Dactolisib_Trametinib_combination dataset
#' 
#' Dataset containing the effect of the combination of Dactolisib and Trametinib on the number of BT-20 cell line cell counts.
#' @name Dactolisib_Trametinib_combination
#' @docType data
#'
#'@usage data(Dactolisib_Trametinib_combination)
#'
#'@details
#'  The data consists of a dataframe named Dactolisib_Trametinib_combination with the following columns:
#'  \itemize{
#'   \item Cell.line: name of the cell line used in the experiment. BT-20 in this example.
#'   \item Small.Molecule.Name: Name of one of the drug used in the experiment. Dactolisib in this case.
#'   \item CONC: concentration of the drug specified in Small.Molecule.Name column.
#'   \item CONC.units: units of drug concentrations in CONC. Numeric column.
#'   \item Small.Molecule.2.Name: Name of the second drug used in the experiment. Trametinib in this case.
#'   \item CONC2: concentration of the drug specified in Small.Molecule.2.Name column. Numeric column.
#'   \item CONC2.units: units of drug concentrations in CONC2.
#'   \item Biological.Replicate: only 1 in this example. Numeric column.
#'   \item Replicate: Technical Replicate. Numeric column.
#'   \item Cell_Count_0: Total cell counts before treatment.
#'   \item Viable.cells: Total cell count after treatment
#'   \item Control: Total control cell count.
#'   \item Normalized_GR: Normalized Growth Rate Inhibition value.
#'   \item Time: Time where the viable cells are counted after treatment.
#'   \item Type: drug sensitive cell type. Number 0 in this case.
#'   }
#'
#'  
#' 
#'@examples
#'\dontrun{
#'  data(Dactolisib_Trametinib_combination)
#'  print(Dactolisib_Trametinib_combination)
#'}
#'@keywords datasets
NULL

#' Dactolisib_Trametinib_rates dataset
#' 
#' Dataset containing the effect of the combination of Dactolisib and Trametinib on the birth rate of drug sensitive BT-20 cell line.
#' @name GD
#' @docType data
#'
#'@usage data(Dactolisib_Trametinib_rates)
#'
#'@details
#'  The data consists of a dataframe named GD with the following columns:
#'  \itemize{
#'   \item Cell.line: name of the cell line used in the experiment. BT-20 in this example.
#'   \item CONC: concentration of the drug 1. Numeric column.
#'   \item CONC2: concentration of the drug 2. Numeric column.
#'   \item Type: drug sensitive cell type. Number 0 in this case.
#'   \item Net_growth: net growth rates of the cells under different drug concentrations.
#'   \item Birth_rate: birth rate of the sensitive cell line influenced by different drug concentrations.
#'   \item Death_rate: Constant death rate.
#' }
#'
#'  
#' 
#'@examples
#'\dontrun{
#'  data(Dactolisib_Trametinib_rates)
#'  head(GD)
#'}
#'@keywords datasets
NULL

#' AU565 dataset
#' 
#' Dataset containing the effect of a DNA-alkylator on the viable cell counts of AU565 cell line.
#' @name AU565_dataset
#' @docType data
#'
#'@usage data(AU565_dataset)
#'
#'@details
#'  The data consists of a dataframe named GD with the following columns:
#'  \itemize{
#'   \item Cell.line: name of the cell line used in the experiment. AU565 in this example.
#'   \item CONC: concentration of the drug 1. Numeric column.
#'   \item Type: drug sensitive cell type. Number 0 in this case.
#'   \item Viable.cells: Total cell count.
#'   \item Replicate: Technical replicate
#' }
#' 
#'@examples
#'\dontrun{
#'  data(AU565_dataset)
#'  head(AU565_dataset)
#'}
#'@keywords datasets
NULL
