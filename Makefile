all:
	Rscript ComputeMFP.R
	@perl -e "print \"Average number of features to detect anomaly: \";"; awk -F',' '{ sum += $$4 } END { if (NR > 0) print sum / (NR-1) }' MFP.seq_dropout.shuttle.allBenchR10.csv
