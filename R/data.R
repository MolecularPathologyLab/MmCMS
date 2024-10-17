#' Mouse CMS template (option A)
#' @details CMS prediction template for NTP method. Marker genes were
#' identified by converting human CMS template from \href{https://www.nature.com/articles/s41598-017-16747-x}{CMScaller} to mouse orthologe using biomaRt.
#' \code{templates$probe} refers to gene symbols.
#' @references Eide PW, Bruun J, Lothe RA, Sveen A. CMScaller:
#' an R package for consensus molecular subtyping of colorectal cancer pre-clinical models. Scientific reports (2017)
#' @references Guinney J, Dienstmann R, Wang X, de Reynies A, Schlicker A, Soneson C, et al. The consensus molecular subtypes of colorectal cancer. Nat Med (2015)
#' @references Hoshida Y. Nearest Template Prediction: A Single-Sample-Based Flexible Class Prediction with Confidence Assessment. PLoS ONE (2010)
#' @references Durinck S, Spellman P, Birney E, Huber W (2009). “Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt.” Nature Protocols, 4, 1184–1191.
#' @examples
#' emat <- TestData_gemm
#' re <- MmCMS(emat, templates=MmCMS::template.CMS.A, Genesets = c("template.CMS.A"), seed=120)
"template.CMS.A"

#' Mouse CMS template (option B)
#' @details CMS prediction template for NTP method. Marker genes were
#' identified via a biological knowledge-based approach instead of having sole emphesis on individual genes.
#' ssGSEA was performed on TCGA RNA-sequencing data from original \href{https://pubmed.ncbi.nlm.nih.gov/26457759/}{CMSclassifier}
#' with assigned CMS subtypes, to provide GO_BP scores for each sample.
#' Then CMS-related GO_BP terms were converted to mouse GO_BP terms orthologe using msigdbr.
#' Instead of individual genes, the GO_BP terms used as a template.
#' \code{templates$probe} refers to GO_BP terms.
#' @references Eide PW, Bruun J, Lothe RA, Sveen A. CMScaller:
#' an R package for consensus molecular subtyping of colorectal cancer pre-clinical models. Scientific reports (2017)
#' @references Guinney J, Dienstmann R, Wang X, de Reynies A, et al. The consensus molecular subtypes of colorectal cancer. Nat Med (2015)
#' @references Hoshida Y. Nearest Template Prediction: A Single-Sample-Based Flexible Class Prediction with Confidence Assessment. PLoS ONE. 2010;5:e15543.
#' @references Liberzon et al. (2015) The Molecular Signatures Database Hallmark Gene Set Collection. DOI:https://doi.org/10.1016/j.cels.2015.12.004
#' @references Hänzelmann S., Castelo R. and Guinney J. GSVA: gene set variation analysis for microarray and RNA-Seq data. BMC Bioinformatics, 14:7, 2013.
#' @examples
#' emat <- TestData_gemm
#' re <- MmCMS(emat, templates=MmCMS::template.CMS.B, Genesets = c("template.CMS.B"), seed=120)
"template.CMS.B"

#' Mouse CMS template (option C)
#' @details CMS prediction template for NTP method. Marker genes were identified through a more complicated approach.
#' A combined pathway system encompassing multiple biological and histological signalling cascades.
#'\code{templates$probe} refers to gene sets.
#' @references Eide PW, Bruun J, Lothe RA, Sveen A. CMScaller:
#' an R package for consensus molecular subtyping of colorectal cancer pre-clinical models. Scientific reports (2017)
#' @references Guinney J, Dienstmann R, Wang X, de Reynies A, et al. The consensus molecular subtypes of colorectal cancer. Nat Med (2015)
#' @references Hoshida Y. Nearest Template Prediction: A Single-Sample-Based Flexible Class Prediction with Confidence Assessment. PLoS ONE. (2010)
#' @references Liberzon et al. The Molecular Signatures Database Hallmark Gene Set Collection. Cell Syst (2015)
#' @references Petitprez, F., Lévy, S., Sun, C.-M., Meylan, M., et al. The murine Microenvironment Cell Population counter method to estimate abundance of tissue-infiltrating immune and stromal cell populations in murine samples using gene expression. Genome Med (2020)
#' @references Hänzelmann S., Castelo R. and Guinney J. GSVA: gene set variation analysis for microarray and RNA-Seq data. BMC Bioinformatics, 14:7, 2013.
#' @examples
#' emat <- TestData_gemm
#' re <- MmCMS(emat, templates=MmCMS::template.CMS.C, Genesets = c("template.CMS.C"), seed=120)
"template.CMS.C"

#' Expression profile of GEMMs
#' @details 18 GEMM data from Cancer Cell. 2019 Sep 16;36(3):319-336.e7.doi: 10.1016/j.ccell.2019.08.003.
"TestData_gemm"

#' GO_BP geneset for mouse
#' @details GO_BP terms from msigdb
"GO.BP.list_mus"

#' option.C geneset for mouse
#' @details gene set includes a combined pathway system
"option.C.geneset.list_mus"
