
#-------------------------------------------------------------------#
#                                                                   #
#   Prepare the Gamma-shape by DAF file using neutral simulations   #
#                                                                   #
#-------------------------------------------------------------------#


#####------------------- CHANGE THESE LINES! --------------------####

MS				=	./ms
SIMUPOP_BACKWARD		=	./backward.py -m Tennessen_CEU
PRESENT_DIPLOID_POPULATION_SIZE	=	512000
#PRESENT_DAF_LIST		=	$(shell seq 0.05 0.01 0.95)
#PRESENT_DAF_LIST		=	$(shell seq 0.05 0.05 0.95)
#the variation in gamma-shape values is higher for the low daf regime.
#also, the higher the sample size the lower the resolution needed
#this should be enough for n=3,000:
PRESENT_DAF_LIST		=	0.05 0.07 0.10 0.20 0.50 0.95
SIM_NUM_REPLICATIONS		=	1000
#sample size is number of haplotypes (so twice the number of individuals)
SAMPLE_SIZE			=	6000

#####------------------------------------------------------------####




MS_DERIVED_SAMPLE_SIZE		=	$(shell perl -e '$$i=int($(SAMPLE_SIZE)*$(DAF));print"$$i";')
MS_DERIVED_INIT_DIPL_POPSIZE	=	$(shell perl -e '$$i=$(PRESENT_DIPLOID_POPULATION_SIZE)*$(DAF);print"$$i";')


sim:
	echo -n > sds_input.gamma_shapes; \
	\
	$(foreach daf, $(PRESENT_DAF_LIST), \
		make my_simulation0 DAF=$(daf); \
	) \


my_simulation0:
	echo -n > tmp_stats.tab; \
	\
	perl -e 'for($$i=1;$$i<=$(SIM_NUM_REPLICATIONS);$$i++){system("make my_simulation1");}'; \
	\
	cut -f2 tmp_stats.tab \
	| perl -e '@r=();while($$l=<STDIN>){chomp($$l);push(@r,$$l);} $$ave=ave(@r); print"$(DAF)\t$$ave\n";sub ave {my(@v)=@_;my($$c)=scalar(@v);my($$t)=0;$$t+=$$_ for @v; return $$c?$$t/$$c:0;}' \
	>> sds_input.gamma_shapes; \
	\
	rm -f tmp_stats.tab; \


my_simulation1:
	make my_simulation_prepare_trajectory; \
	make my_simulation_trajectory_2_trees; \
	make my_simulation_trees_2_tip_branch_stats; \
	make my_simulation_clean; \


my_simulation_prepare_trajectory:
	$(SIMUPOP_BACKWARD) -te 0 -ts 100000 -ps 0 -g 100000 -f $(DAF) \
	| perl -e '$$l=<STDIN>;chomp($$l);@r=split(/\t/,$$l);$$t0=0.0;$$t1=$$r[0];$$n=$$r[1];$$freq=$$r[5];$$sigma=$$r[3]-1.0;print"$$t0\t$$t1\t$$n\t$$freq\t$$sigma\n";while($$l=<STDIN>){chomp($$l);@r=split(/\t/,$$l);$$t0=$$t1;$$t1=$$r[0];$$n=$$r[1];$$freq=$$r[5];$$sigma=$$r[3]-1.0;print"$$t0\t$$t1\t$$n\t$$freq\t$$sigma\n";}' \
	| perl -e '$$l=<STDIN>;chomp($$l);@r=split(/\t/,$$l);$$n0=$$r[2];$$tscale=4.0*$$n0;$$t0=$$r[0]/$$tscale;$$t1=$$r[1]/$$tscale;$$n=$$r[2]/$$n0;$$freq=$$r[3];$$sigma=$$r[4];print"$$t0\t$$t1\t$$n\t$$freq\t$$sigma\n";while($$l=<STDIN>){chomp($$l);@r=split(/\t/,$$l);$$t0=$$r[0]/$$tscale;$$t1=$$r[1]/$$tscale;$$n=$$r[2]/$$n0;$$freq=$$r[3];$$sigma=$$r[4];print"$$t0\t$$t1\t$$n\t$$freq\t$$sigma\n";}' \
	| perl -e '$$l=<STDIN>;chomp($$l);@r=split(/\t/,$$l);$$t0=$$r[0];$$t1=$$r[1];$$n=$$r[2];$$p=$$r[3];$$sigma=$$r[4];while($$l1=<STDIN>){chomp($$l1);@r=split(/\t/,$$l1);if($$n!=$$r[2]||$$p!=$$r[3]){print"$$t0\t$$t1\t$$n\t$$p\t$$sigma\n";$$l=$$l1;$$t0=$$r[0];$$t1=$$r[1];$$n=$$r[2];$$p=$$r[3];$$sigma=$$r[4];} else{$$t1=$$r[1];}} print"$$t0\t$$t1\t$$n\t$$p\t$$sigma\n"; print"$$t1\t999.0\t$$n\t0\t0\n";' \
	| cut -f1-4 \
	> tmp_sim_daf_traj.tab; \


my_simulation_trajectory_2_trees:
	cat tmp_sim_daf_traj.tab \
	| perl -e '$$l=<STDIN>;chomp($$l);@r=split(/\t/,$$l);$$AF0=$$r[3];while($$l=<STDIN>){chomp($$l);@r=split(/\t/,$$l);$$Ta=$$r[0]/$$AF0;$$Na=$$r[2]*$$r[3]/$$AF0;$$Tb=$$r[0]/(1-$$AF0);$$Nb=$$r[2]*(1-$$r[3])/(1-$$AF0);print"-eN\t$$Ta\t$$Na\t$$Tb\t$$Nb\t$$r[0]\n";}' \
	| cut -f 1,2,3 \
	> tmp_sim_msargs.d; \
	\
	$(MS) $(MS_DERIVED_SAMPLE_SIZE) 1 -p 20 -T -f tmp_sim_msargs.d \
	| tee /dev/null \
	| tail -n1 | sed 's/\([0-9][0-9]*\)\:/D&/g' \
	> tmp_sim_newick_tree.d; \


my_simulation_trees_2_tip_branch_stats:
	cat tmp_sim_newick_tree.d \
	| perl -e 'use List::Util qw(shuffle); $$lin0=<STDIN>;chomp($$lin0);$$k=1; @r=(1..$(MS_DERIVED_SAMPLE_SIZE));@rs=shuffle(@r);$$list_length=@r;$$i0=0;while($$k>1){$$i=int($$list_length/$$k);$$lin=$$lin0;for($$j=$$i0;$$j<$$i0+$$i;$$j++){$$lin=~s/D$$rs[$$j]\:/M$$rs[$$j]\:/;} my_func($$lin,$$k);$$i0+=$$i;--$$k;$$list_length-=$$i; } $$lin=$$lin0; for($$j=$$i0;$$j<@rs;$$j++){$$lin=~s/D$$rs[$$j]\:/M$$rs[$$j]\:/;} my_func($$lin,$$k); sub my_func {($$lin,$$kk)=@_;$$l=$$lin;while(1){if($$l=~m/(.*)\(([M|Y|A|D|X][0-9]*)\:(\d*.?\d*)\,([M|Y|A|D|X][0-9]*)\:(\d*.?\d*)\)(.*)/){$$n1=$$2;$$n2=$$4;$$e1=$$3;$$e2=$$5;$$l1=$$1;$$l2=$$6;$$d1=($$n1=~m/[XY](\d*)/)?$$1:1;$$d2=($$n2=~m/[XY](\d*)/)?$$1:1;$$d=$$d1+$$d2; $$isin1=($$n1=~m/[MY]/); $$isin2=($$n2=~m/[MY]/);$$isin=($$isin1 and $$isin2)?1:0;print"$$d1\t$$e1\t$$n1\t$$kk\n$$d2\t$$e2\t$$n2\t$$kk\n"; if($$isin){$$l=$$l1."Y".$$d.$$l2;} elsif($$isin1){$$l=$$l1.$$n1.$$l2;} elsif($$isin2){$$l=$$l1.$$n2.$$l2;} else{$$l=$$l1."X"."0".$$l2;}} else{last;}} return;}' \
	| perl -e 'while($$l=<STDIN>){chomp($$l);@r=split(/\t/,$$l);if($$l==1){print"$$l\n";}}' | grep M | cut -f 2 \
	| perl -e 'while($$l=<STDIN>){chomp($$l);$$v=4.0*$(MS_DERIVED_INIT_DIPL_POPSIZE)*$$l;print"$$v\n";}' \
	| perl -e '@r=();while($$l=<STDIN>){chomp($$l);push(@r,$$l);} $$ave=ave(@r); $$sd=std($$ave,@r); $$shape=($$ave*$$ave)/($$sd*$$sd); print"$(DAF)\t$$shape\n";sub ave {my(@v)=@_;my($$c)=scalar(@v);my($$t)=0;$$t+=$$_ for @v; return $$c?$$t/$$c:0;} sub std {my($$ave,@v)=@_;my($$c)=scalar(@v);my($$s)=0;$$s+=($$_-$$ave)**2 for @v; return $$c?sqrt($$s/$$c):0;}' \
	>> tmp_stats.tab; \


my_simulation_clean:
	rm -f 	tmp_sim_daf_traj0.tab \
		tmp_sim_daf_traj.tab \
		tmp_sim_msargs.d \
		tmp_sim_newick_tree.d \
	; \




#----------------------------------------#
#                                        #
#   Run SDS on the example input files   #
#                                        #
#----------------------------------------#

# This example is based on the simulation presented in Fig1 in our paper.
#  (selection from standing variation starting 100 generations ago; 
#   Tennessen CEU demographic model; present DAF=0.7; s=0.05; n=3,000 individuals)
# 
# Note that the singletons file used here shows only part of the total singletons in the region,
# but there is no reason to do it. In practice I used a file for all singletons in a chromosome.
# I then could break the chromosome-wide test-SNP file into small files and run in parallel
# (using the same chromosome-wide sinlgetons file).
#
# If you uncomment the "DEBUG_MODE=TRUE" line in compute_SDS.R you'll get a couple of figures and more stats
# (the program will halt after each figure -- press any key to continue)

run_example:
	./compute_SDS.R \
		example.singletons \
		example.testsnp \
		example.observability \
		example.boundaries \
		example.gamma_shapes \
		1e-6 \
		10 \
	> example.output; \

