#! /bin/bash

for t in 1 2 3 4; do
    for a in 4 6 5 7; do 
	for n in `seq $a 4 36`; do 
	    go test -v -run=SetExternalParam -bench=2N.* -cutoff=$n -blocksize=$n -benchtime=4s >> 2N.out
	    awk 'BEGIN{print "# n, co/bs Rec. B.U."}
             /^Benchmark2NStable/{s=$3; next}
             /^Benchmark2NBUStable/{
                 gsub(".*Stable","",$1);
                 gsub("K","000",$1);
                 printf "%s %i %.2f %.2f\n", $1, co, $3/1000, s/1000; next
             }
             /Cutoff/{co=$4; bs=$8; next;}' 2N.out | sed -e "s/,/./g" > 2N
	    gnuplot measure.gp
	done
    done

    for a in 4 6 5 7; do 
	for n in `seq $a 4 36`; do 
	    go test -v -run=SetExternalParam -bench=N2.* -cutoff=$n -blocksize=$n -benchtime=4s >> N2.out
	    awk 'BEGIN{print "# n, co/bs Rec. B.U."}
             /^BenchmarkN2Stable/{s=$3; next}
             /^BenchmarkN2BUStable/{
                 gsub(".*Stable","",$1);
                 gsub("K","000",$1);
                 printf "%s %i %.2f %.2f\n", $1, co, $3/1000, s/1000; next
             }
             /Cutoff/{co=$4; bs=$8; next;}' N2.out | sed -e "s/,/./g" > N2
	    gnuplot measure.gp
	done
    done
done

