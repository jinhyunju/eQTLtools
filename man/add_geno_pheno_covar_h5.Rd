% Generated by roxygen2 (4.0.2): do not edit by hand
\name{add_geno_pheno_covar_h5}
\alias{add_geno_pheno_covar_h5}
\title{Adding genotypes and phenotypes to hdf5 file}
\usage{
add_geno_pheno_covar_h5(file_name, phenotypes = NULL, genotypes = NULL,
  covars = NULL)
}
\arguments{
\item{file_name}{Name of the HDF5 file to add the information to}

\item{phenotypes}{Phenotype matrix to save}

\item{genotypes}{Genotype matrix to save}

\item{covars}{Covariate matrix to save}
}
\value{
Saves the given objects into the hdf5 generated by the create_eqtl_input_h5() function.
}
\description{
Use in-memory objects to create genotypes and phenotypes datasets in a given hdf5 file
}
\keyword{HDF5}
\keyword{eQTL,}

