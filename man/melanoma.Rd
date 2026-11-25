\name{melanoma}
\alias{melanoma}
\docType{data}
\title{Malignant melanoma data prepared for competing-risks model comparison}
\description{
Analysis-ready data for the 205 patients in the melanoma application.
The original measurements are from patients undergoing surgery for malignant
melanoma at the University Hospital of Odense, Denmark, during 1962--1977.
}
\usage{melanoma}
\format{
A data frame with 205 rows and 7 variables:
\describe{
\item{id}{Sequential row identifier in the original dataset, from 1 to 205.}
\item{time}{Observed follow-up in years, obtained by dividing days by 365.25.}
\item{status}{Event code: 0, censored (134); 1, melanoma death (57);
2, death from other causes (14).}
\item{age}{Age at surgery, in years.}
\item{sex}{Factor with levels \code{female} and \code{male}, in that order.}
\item{thick}{Natural logarithm of one plus tumour thickness in millimetres.}
\item{ulcer}{Factor with levels \code{absent} and \code{present}, in that order.}
}
}
\details{
This dataset is processed from \code{boot::melanoma}, version 1.3-32.
All 205 original rows are retained. Original status values 1, 2 and 3 are
mapped to 1, 0 and 2, respectively. The original calendar year is omitted.
The time conversion, transformation and factor reference levels match the
paper's analysis. No patient information was inferred or imputed.
This is distinct from the unrelated dataset of the same name in lattice.
}
\source{
The \pkg{boot} package, distributed with License: Unlimited.
\url{https://stat.ethz.ch/R-manual/R-devel/library/boot/html/melanoma.html}
}
\references{
Andersen, P. K., Borgan, O., Gill, R. D. and Keiding, N. (1993).
\emph{Statistical Models Based on Counting Processes}. Springer-Verlag.
}
\examples{
data("melanoma", package = "csplcsh")
table(melanoma$status)
head(melanoma)
}
\keyword{datasets}
