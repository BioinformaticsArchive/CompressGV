#!/bin/bash

set -e
TESTMODE=0;

##########
# Get test data
##########

echo -e "----------\nCreating test data - a large amount of data needs to be downloaded if not already present (>2.5G) \n----------";

# Check that the correct HumVar dataset is available
cd ./data;

function downloadAndCheck { #$1 = local file prefix; $2 = [gz|bz2]; $3 = URL; $4 = error message
	([ -f $1.tar.$2 ] && md5sum -c $1.md5 >/dev/null 2>&1) \
		|| ( \
			echo "Downloading $3" && wget --progress=bar:force -O $1.tar.$2 $3 2>&1 | tail -f -n +12 \
			&& ( \
				md5sum -c $1.md5 >/dev/null 2>&1 || (echo "ERROR: $4" && exit 1)
			) \
		);
}

downloadAndCheck humvar gz "ftp://genetics.bwh.harvard.edu/pph2/training/humvar-2011_12.predictions.tar.gz" "HumVar dataset does not match expected hash";

echo -n "HumVar dataset exists locally and matches expected hash; extracting... ";
mkdir -p humvar;
cd ./humvar;
tar -xzf ../humvar.tar.gz humvar-2011_12.deleterious.humvar.output humvar-2011_12.neutral.humvar.output;
echo "done.";

echo -n "Finding proteins with at least 3 variants in each of deleterious and neutral sets... ";

function atLeast3 { #$1 = [deleterious|neutral]
	tail -n +2 humvar-2011_12.$1.humvar.output | awk '{print $1}' | sort | uniq -c | awk '$1>2 {print $2}'
}

comm -12 <(atLeast3 deleterious) <(atLeast3 neutral) > toTest.log;
echo `< toTest.log wc -l` "found.";

echo "Extracting variants for included proteins... ";
mkdir -p deleterious;
mkdir -p neutral;

function extractVariants { #$1 = [deleterious|neutral]; $2 = variant accession
	tail -n +2 humvar-2011_12.$1.humvar.output | awk -v variant="$2" '$1==variant {print $3$2$4}' > $1/$2;
}

n=`< toTest.log wc -l`;
i=0;
if [ $TESTMODE -eq 1 ]
then
	echo "[TEST MODE] Skipping extraction and keeping existing files";
else
	for v in $(cat toTest.log)
	do
		echo -ne "$i of $n\r";
		extractVariants deleterious $v;
		extractVariants neutral $v;
		i=`expr $i + 1`;
	done
	echo -e "Extracted.\r";
fi

cd ..;
echo -n "Checking for precomputed alignments... ";
if [ $TESTMODE -eq 1 ]
then
	echo "[TEST MODE] Assuming present and valid.";
else
	downloadAndCheck alignments bz2 "ftp://genetics.bwh.harvard.edu/pph2/bundled/polyphen-2.2.2-alignments-mlc-2011_12.tar.bz2" "Precomputed alignments do not match expected hash";
	echo "done.";
fi

if [ -d alignments/polyphen-2.2.2/precomputed/alignments ]
then
	echo "Alignments directory exists; skipping extraction.";
else
	mkdir -p alignments;
	cd alignments;
	cat ../humvar/toTest.log | xargs -i echo "polyphen-2.2.2/precomputed/alignments/{}.aln" > toExtract.log;
	echo "Precomputed alignments exist locally and match expected hash; extracting... ";
	
	tar -xjvf ../alignments.tar.bz2 --files-from toExtract.log > extracted.log &
	pid=$!;

	n=`< toExtract.log wc -l`;
	i=0
	while [ $i -lt $n ]
	do
		echo -ne "$i of $n\r";
		i=`< extracted.log wc -l`;
		sleep 2;
	done
	
	kill $pid;
	
	echo -e "Extracted.\r";
	cd ..;
fi

cd alignments;
mkdir -p fasta;
cd fasta;
echo -n "Downloading FASTA sequences from UniProt... ";
cat ../../humvar/toTest.log | xargs -n 1 -i -P 20 bash -c "[ -f {}.fasta ] || wget -q http://www.uniprot.org/uniprot/{}.fasta";
echo "done.";
cd ..;

echo "Combining FASTA sequences with PolyPhen precomputed alignments as they do not include the human protein...";
mkdir -p combined;
cd combined;
n=`ls ../fasta | fgrep .fasta | wc -l`;
i=0;
for g in $(ls ../fasta | fgrep .fasta | cut -b 1-6)
do
	echo -ne "$i of $n\r";
	if [ ! -e ${g}.aln ]
	then
		# Extract the alignments from polyphen-2.2.2 precomputed into unlabeled FASTA
		tail -n +3 ../polyphen-2.2.2/precomputed/alignments/${g}.aln | fgrep -v UPI | awk '{print ">\n"$NF}' > ${g}.aln.fa;
		# Add the human FASTA
		clustalo --p1 ../fasta/${g}.fasta --p2 ${g}.aln.fa > ${g}.aln.hum.fa
		# Remove the potential error message caused by on leading and/or trailing residues in the alignment
		head -n 1 ${g}.aln.hum.fa | grep -P '^Transfer:' >/dev/null && tail -n +2 ${g}.aln.hum.fa > ${g}.tmp && mv ${g}.tmp ${g}.aln.hum.fa;
		# Back to alignments only
		../../../../fastaToSingleLine.sh ${g}.aln.hum.fa > ${g}.aln
		# Clean up
		rm ${g}.aln.fa ${g}.aln.hum.fa
	fi
	i=`expr $i + 1`;
done
echo -e "Done.          \r";
cd ../../../;

echo -e "----------\nTest data created\n----------";
echo -e "----------\nRunning cross-validation\n----------";

if [ ! -f ../bin/grantham ]
then
	echo "Compiled binary does not exist; attempting to create before running analysis.";
	currDir=`pwd`;
	cd ..;
	make;
	cd $currDir;
	echo "Binary compiled; commencing validation.";
fi

mkdir -p data/outcome;
proc=`< /proc/cpuinfo fgrep processor | wc -l | xargs -i expr 2 \* {}`;
echo "Running ${proc} processes simultaneously";
cat data/humvar/toTest.log | xargs -n 1 -P $proc ./validate-one.sh > /dev/null &

n=`< data/humvar/toTest.log wc -l`;
i=0
e=0
while [ $i -lt $n ] 
do
	i=`ls data/outcome | fgrep verbose | wc -l`
	e=`ls data/outcome | fgrep verbose.err | wc -l`
	echo -ne "$i of $n (${e} errors)\r";
done

echo -e "Validation complete; analysing...\r";

if [ $e -gt 0 ]
then
	echo "Errors when validating; check data/outcome/*.err";
	exit 1;
fi

mkdir -p data/analysis
cd data/outcome;
echo -e 'TP\nTN\nFP\nFN' | xargs -n 1 -i bash -c "ls | fgrep verbose | xargs fgrep --text {} | awk '{print \$NF}' > ../analysis/{}";
cd ../analysis;
cat TP FN > P;
cat TN FP > N;
cd ../../;

R --vanilla < ./performace.R 2>/dev/null | tail -n 7
