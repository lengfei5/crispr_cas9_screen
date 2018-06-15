### script to run Georg's screen analysis in cluster
nb_cores=20;
qsub -q public.q -o ${PWD}/logs -j yes -pe smp $nb_cores -cwd -b y -shell y -N screen "bash run_all.sh;"
