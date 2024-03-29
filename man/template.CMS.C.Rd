% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{template.CMS.C}
\alias{template.CMS.C}
\title{Mouse CMS template (option C)}
\format{
An object of class \code{data.frame} with 77 rows and 2 columns.
}
\usage{
template.CMS.C
}
\description{
Mouse CMS template (option C)
}
\details{
CMS prediction template for NTP method. Marker genes were identified through a more complicated approach.
A combined pathway system encompassing multiple biological and histological signalling cascades.
\code{templates$probe} refers to gene sets.
}
\examples{
emat <- TestData_gemm
re <- MmCMS(emat, templates=MmCMS::template.CMS.C, Genesets = c("template.CMS.C"), seed=120)
}
\references{
Eide PW, Bruun J, Lothe RA, Sveen A. CMScaller:
an R package for consensus molecular subtyping of colorectal cancer pre-clinical models. Scientific reports (2017)

Guinney J, Dienstmann R, Wang X, de Reynies A, et al. The consensus molecular subtypes of colorectal cancer. Nat Med (2015)

Hoshida Y. Nearest Template Prediction: A Single-Sample-Based Flexible Class Prediction with Confidence Assessment. PLoS ONE. (2010)

Liberzon et al. The Molecular Signatures Database Hallmark Gene Set Collection. Cell Syst (2015)

Petitprez, F., Lévy, S., Sun, C.-M., Meylan, M., et al. The murine Microenvironment Cell Population counter method to estimate abundance of tissue-infiltrating immune and stromal cell populations in murine samples using gene expression. Genome Med (2020)

Hänzelmann S., Castelo R. and Guinney J. GSVA: gene set variation analysis for microarray and RNA-Seq data. BMC Bioinformatics, 14:7, 2013.
}
\keyword{datasets}
