 # for n in `ls F32A_1/*.seq F32A_2/*.seq`; do echo '>'$n; cat $n; done > F32.fna
 # for n in `ls MutL_1/*.seq MutL_2/*.seq`; do echo '>'$n; cat $n; done > MutL.fna
 # for n in `ls MutSdel_1/*.seq MutSdel_2/*.seq`; do echo '>'$n; cat $n; done > MutSdel.fna
 # for n in `ls Parent_1/*.seq Parent_2/*.seq`; do echo '>'$n; cat $n; done > Parent.fna



./SilentDNA.pl F32.fna > F32.csv &
./SilentDNA.pl MutSdel.fna > MutSdel.csv &
./SilentDNA.pl MutL.fna > MutL.csv &
./SilentDNA.pl Parent.fna > Parent.csv &
