#!/bin/sh

if [ $# -ne 8 ]; then
    if [ $# -ne 9 ]; then
	echo 'usage:'
	echo '   ./slide [-G|-H] [data file] [#control] [#case] [w] [tmpdir] [#sampling] [seed] [mapfile(option)]\n'
	exit 0
    fi
fi

mkdir -p $6
./slide_1prep $1 $2 $3 $4 $5 $6/slide
./slide_2run $6/slide $6/slide.maxstat $7 $8
./slide_3sort $6/slide.sorted $6/slide.maxstat
./slide_4correct -p $6/slide.sorted $6/slide.txt_pointwisep $6/slide.pvalues $9
