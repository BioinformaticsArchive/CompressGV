if [ ! -f ../Debug/grantham ]
then
	echo "Compiled binary does not exist; attempting to create before running demo.";
	currDir=`pwd`;
	cd ..;
	make;
	cd $currDir;
	echo "Binary compiled; running demo.";
fi

../Debug/grantham ./aln ./del ./neut ./novel
