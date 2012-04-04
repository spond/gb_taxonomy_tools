These are four simple utilities which perform the following manipulations and visualization tasks on GenBank 
taxonomic information.

+ `gid-taxid` : convert a list of GenBank IDs and associated counts into the list of tripets: genbank id, taxonomy id, count. 
It requires access to (quite large) mapping files maintained by GenBank, e.g. ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.zip
which are tab separated lists of  `gid taxid count`, e.g. the input line `160338813  160` is output as `160338813  436308	160`
Try running it on as `$gid-taxid tests/data/test.gid path/to/gi_taxid_nucl.dmp` The result should be as in tests/data/test.taxid 

+ `taxonomy-reader`