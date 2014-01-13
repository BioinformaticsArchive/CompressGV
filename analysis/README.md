./run.sh will download all necessary test data, extract relevant files, and build the required alignments in order to test CompressGV. Cross-validation will be performed and results analysed with R. Tested in Ubuntu 13.10.

# Requirements

- Internet connection with >2.5G download limit remaining
- wget
- md5sum
- clustalo
- R + ROCR library

## Notes

- The file *humvar_errors* contains a list of proteins / variants that appear to be errors in the HumVar dataset as the wild-type in the variant does not match that in the Uniprot sequence. Cursory analysis appears to show inversion of the wild-type / variant amino acids in the HumVar data