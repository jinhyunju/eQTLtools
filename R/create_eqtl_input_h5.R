#' Creating eQTL input hdf5 file templates.
#'
#' Generates an hdf5 file with the basic groups for eQTL analysis.
#'
#' @param file_name The name of the file to be created
#'
#' @return HDF5 file with groups phenotypes, genotypes, and covars with subgroups col_info and row_info will be created in the current directory.
#' @keywords eQTL, HDF5
#'
#' @import rhdf5
#' @export
create_eqtl_input_h5 <- function(file_name){
	level1.groups <- c("phenotypes", "genotypes", "covars")
    
	h5createFile(file_name)
	
    for(l1 in 1:length(level1.groups)){
        h5createGroup(file_name, level1.groups[l1])
        h5createGroup(file_name, paste(level1.groups[l1], "col_info", sep = "/"))
        h5createGroup(file_name, paste(level1.groups[l1], "row_info", sep = "/"))
    }
	
}