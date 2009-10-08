cnF2freq <-
function(ped.file, chr.file, os = "unix", compile = FALSE) {
## --------------------------------------------- ##
##        Run Compiled cnF2freq Software         ##
##      Xia.Shen@lcb.uu.se ---2009-09-14---      ##
## --------------------------------------------- ##
	if (compile) {
		## Compile cnF2freq
		system(paste("g++ -I . -O2 cnF2freq.cpp -o cnF2freq_", os, sep = ""))
	}
	## Rename input files for cnF2freq
	copy <- file.copy(ped.file, "pedi_temp_xiashen.ric")
	copy <- file.copy(chr.file, "chro_temp_xiashen.ric")
	rename <- file.rename("pedi_temp_xiashen.ric", "ped.ric")
	rename <- file.rename("chro_temp_xiashen.ric", "chr16.ric")
	## Run cnF2freq
	if (os == "windows") {
        system(paste("cnF2freq_", os, sep = ""))    
    } else {
        system(paste("./cnF2freq_", os, sep = ""))
    }
##END
}

