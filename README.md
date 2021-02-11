

scanMatchedFilter.py scans the whole genome with matched filter and applies the SVM model using training data provided with the code to predict enhancers and promoters in a cell-type specific manner.

python scanMatchedFilter.py <fileList> <metaProfileList> <chrNameList> <peakFileList> <opPrefix>

where:

	<fileList> is a file with the list of chromatin signals in the format (2 column tab delimited with experimental dataset name in column 1 and filename in column 2). The chromatin signals are in log foldchange signal signal enrichment over control.

	<metaProfileList> is the list with training profiles. This is a tab-delimited 2 column file with the first column containing experimental dataset name (for example, H3K4me1) and the 2nd column containing file name with metaprofile (provided in subdirectory metaprofiles).

	<chrNameList> is the list with chromosome names. The first column contains chromosome name and 2nd column contains length of chromosome.
	
	<peakFileList> is the file with chromatin peaks. These regions are removed during background model fitting.
	
	<positiveScores> is the file containing the scores for all training positives, provided in the training data directory.
	
	<negativeScores> is the file containing the scores for all training negatives, provided in the training data directory.

	<opPrefix> is the prefix for all output files.

	The final output file test_SVMpredScores.dat contains the SVM scores and predictions.
	
