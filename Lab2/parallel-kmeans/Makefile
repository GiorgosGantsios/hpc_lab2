#  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#    File:         Makefile                                                  */
#    Description:  Makefile for programs running a simple k-means clustering */
#                  algorithm                                                 */
#                                                                            */
#    Author:  Wei-keng Liao                                                  */
#             ECE Department Northwestern University                         */
#             email: wkliao@ece.northwestern.edu                             */
#                                                                            */
#    Copyright (C) 2005, Northwestern University                             */
#    See COPYRIGHT notice in top-level directory.                            */
#                                                                            */
#  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

.KEEP_STATE:

all: seq 


CC          = icx

INCFLAGS    = -I.
OPTFLAGS    = -O4 -DNDEBUG -qopenmp -g
LDFLAGS     = $(OPTFLAGS)


CFLAGS      = $(OPTFLAGS) $(DFLAGS) $(INCFLAGS)


H_FILES     = kmeans.h

COMM_SRC = file_io.c util.c


SEQ_SRC     = seq_main.c   \
              seq_kmeans.c \
	      wtime.c

SEQ_OBJ     = $(SEQ_SRC:%.c=%.o) $(COMM_SRC:%.c=%.o)

$(SEQ_OBJ): $(H_FILES)

seq_main.o: seq_main.c $(H_FILES)
	$(CC) $(CFLAGS) -c $*.c

seq_kmeans.o: seq_kmeans.c $(H_FILES)
	$(CC) $(CFLAGS) -c $*.c

wtime.o: wtime.c
	$(CC) $(CFLAGS) -c $*.c

seq: seq_main
seq_main: $(SEQ_OBJ) $(H_FILES)
	$(CC) $(LDFLAGS) -o $@ $(SEQ_OBJ) $(LIBS)

IMAGE_FILES =   color100.txt   color17695.bin   color17695.nc \
                 edge100.txt    edge17695.bin    edge17695.nc \
              texture100.txt texture17695.bin texture17695.nc

INPUTS = $(IMAGE_FILES:%=Image_data/%)

PACKING_LIST = $(COMM_SRC) $(SEQ_SRC) $(H_FILES) \
               Makefile README COPYRIGHT

dist:
	dist_dir=parallel-kmeans \
	&& rm -rf $$dist_dir $$dist_dir.tar.gz\
	&& mkdir -p $$dist_dir/Image_data \
	&& cp $(PACKING_LIST) $$dist_dir \
	&& cp $(INPUTS) $$dist_dir/Image_data \
	&& tar -cf - $$dist_dir | gzip > $$dist_dir.tar.gz \
	&& rm -rf $$dist_dir

clean:
	rm -rf *.o seq_main \
		core* .make.state              \
		*.cluster_centres *.membership \
		*.cluster_centres.nc *.membership.nc \
		Image_data/*.cluster_centres Image_data/*.membership \
		Image_data/*.cluster_centres.nc Image_data/*.membership.nc

check: all
	# sequential K-means ---------------------------------------------------
	./seq_main -q -b -n 4 -i Image_data/color17695.bin
	./seq_main -q    -n 4 -i Image_data/color100.txt

# New rule to run the executable 12 times and capture time and PSNR
# make run_experiment prefix=test
run_experiment: $(EXECUTABLES)
#	@if [ -z "$(prefix)" ]; then \
#		echo "Error: Please provide a 'prefix' argument for the CSV files."; \
#		exit 1; \
#	fi
#	ORIGINAL
	@rm -f $(prefix)original.csv # Clean up previous run CSVs
	@echo "Total time,PSNR" > $(prefix)original.csv # CSV header
	@total_time_sum=0; total_time_sq_sum=0; count=0; \
	for i in `seq 1 12`; do \
		output=$$(./sobel_orig); \
		time=$$(echo "$$output" | grep "Total time" | awk '{print $$4}'); \
		psnr=$$(echo "$$output" | grep "PSNR" | awk '{print $$9}'); \
		if [ -n "$$time" ]; then \
			echo "$$time,$$psnr" >> $(prefix)original.csv; \
		else \
			echo "Error: could not extract time for execution $$i" >> $(prefix)original.csv; \
		fi; \
	done; \
	python3 calculate_stats.py $(prefix)original
	@echo "Run completed. Results saved in $(prefix)original.csv"
#	################################################################################################
#	LOOP INTERCHANGE
	@echo "Total time,PSNR" > $(prefix)loop_interchange.csv # CSV header
	@total_time_sum=0; total_time_sq_sum=0; count=0; \
	for i in `seq 1 12`; do \
		output=$$(./sobel_Loop_Interchange); \
		time=$$(echo "$$output" | grep "Total time" | awk '{print $$4}'); \
		psnr=$$(echo "$$output" | grep "PSNR" | awk '{print $$9}'); \
		if [ -n "$$time" ]; then \
			echo "$$time,$$psnr" >> $(prefix)loop_interchange.csv; \
		else \
			echo "Error: could not extract time for execution $$i" >> $(prefix)loop_interchange.csv; \
		fi; \
	done; \
	python3 calculate_stats.py $(prefix)loop_interchange
	@echo "Run completed. Results saved in $(prefix)loop_interchange.csv"
#	################################################################################################
# 	LOOP UNROLLING
	@echo "Total time,PSNR" > $(prefix)Loop_Unrolling.csv # CSV header
	@total_time_sum=0; total_time_sq_sum=0; count=0; \
	for i in `seq 1 12`; do \
		output=$$(./sobel_Loop_Unrolling); \
		time=$$(echo "$$output" | grep "Total time" | awk '{print $$4}'); \
		psnr=$$(echo "$$output" | grep "PSNR" | awk '{print $$9}'); \
		if [ -n "$$time" ]; then \
			echo "$$time,$$psnr" >> $(prefix)Loop_Unrolling.csv; \
		else \
			echo "Error: could not extract time for execution $$i" >> $(prefix)Loop_Unrolling.csv; \
		fi; \
	done; \
	python3 calculate_stats.py $(prefix)Loop_Unrolling
	@echo "Run completed. Results saved in $(prefix)Loop_Unrolling.csv"
#	################################################################################################
#	LOOP FUSION
	@echo "Total time,PSNR" > $(prefix)loop_fusion.csv # CSV header
	@total_time_sum=0; total_time_sq_sum=0; count=0; \
	for i in `seq 1 12`; do \
		output=$$(./sobel_Loop_Fusion); \
		time=$$(echo "$$output" | grep "Total time" | awk '{print $$4}'); \
		psnr=$$(echo "$$output" | grep "PSNR" | awk '{print $$9}'); \
		if [ -n "$$time" ]; then \
			echo "$$time,$$psnr" >> $(prefix)loop_fusion.csv; \
		else \
			echo "Error: could not extract time for execution $$i" >> $(prefix)loop_fusion.csv; \
		fi; \
	done; \
	python3 calculate_stats.py $(prefix)loop_fusion
	@echo "Run completed. Results saved in $(prefix)loop_fusion.csv"
#	################################################################################################
#	FUNCTION INLINING
	@echo "Total time,PSNR" > $(prefix)function_inlining.csv # CSV header
	@total_time_sum=0; total_time_sq_sum=0; count=0; \
	for i in `seq 1 12`; do \
		output=$$(./sobel_Function_Inlining); \
		time=$$(echo "$$output" | grep "Total time" | awk '{print $$4}'); \
		psnr=$$(echo "$$output" | grep "PSNR" | awk '{print $$9}'); \
		if [ -n "$$time" ]; then \
			echo "$$time,$$psnr" >> $(prefix)function_inlining.csv; \
		else \
			echo "Error: could not extract time for execution $$i" >> $(prefix)function_inlining.csv; \
		fi; \
	done; \
	python3 calculate_stats.py $(prefix)function_inlining
	@echo "Run completed. Results saved in $(prefix)function_inlining.csv"
#	################################################################################################
#	LOOP INVARIANT CODE MOTION
	@echo "Total time,PSNR" > $(prefix)loop_Invvariant_code_motion.csv # CSV header
	@total_time_sum=0; total_time_sq_sum=0; count=0; \
	for i in `seq 1 12`; do \
		output=$$(./sobel_Loop_Invariant_code_motion); \
		time=$$(echo "$$output" | grep "Total time" | awk '{print $$4}'); \
		psnr=$$(echo "$$output" | grep "PSNR" | awk '{print $$9}'); \
		if [ -n "$$time" ]; then \
			echo "$$time,$$psnr" >> $(prefix)loop_Invvariant_code_motion.csv; \
		else \
			echo "Error: could not extract time for execution $$i" >> $(prefix)loop_Invvariant_code_motion.csv; \
		fi; \
	done; \
	python3 calculate_stats.py $(prefix)loop_Invvariant_code_motion
	@echo "Run completed. Results saved in $(prefix)loop_Invvariant_code_motion.csv"
#	################################################################################################
#	COMMON SUBEXPRESSION ELIMINATION
	@echo "Total time,PSNR" > $(prefix)common_subexpression_elimination.csv # CSV header
	@total_time_sum=0; total_time_sq_sum=0; count=0; \
	for i in `seq 1 12`; do \
		output=$$(./sobel_Common_Subexpression_Elimination); \
		time=$$(echo "$$output" | grep "Total time" | awk '{print $$4}'); \
		psnr=$$(echo "$$output" | grep "PSNR" | awk '{print $$9}'); \
		if [ -n "$$time" ]; then \
			echo "$$time,$$psnr" >> $(prefix)common_subexpression_elimination.csv; \
		else \
			echo "Error: could not extract time for execution $$i" >> $(prefix)common_subexpression_elimination.csv; \
		fi; \
	done; \
	python3 calculate_stats.py $(prefix)common_subexpression_elimination
	@echo "Run completed. Results saved in $(prefix)common_subexpression_elimination.csv"
#	################################################################################################
#	STRENGTH REDUCTION
	@echo "Total time,PSNR" > $(prefix)Strength_Reduction.csv # CSV header
	@total_time_sum=0; total_time_sq_sum=0; count=0; \
	for i in `seq 1 12`; do \
		output=$$(./sobel_Strength_Reduction); \
		time=$$(echo "$$output" | grep "Total time" | awk '{print $$4}'); \
		psnr=$$(echo "$$output" | grep "PSNR" | awk '{print $$9}'); \
		if [ -n "$$time" ]; then \
			echo "$$time,$$psnr" >> $(prefix)Strength_Reduction.csv; \
		else \
			echo "Error: could not extract time for execution $$i" >> $(prefix)Strength_Reduction.csv; \
		fi; \
	done; \
	python3 calculate_stats.py $(prefix)Strength_Reduction
	@echo "Run completed. Results saved in $(prefix)Strength_Reduction.csv"

run_spec: $(EXECUTABLES)
	@echo "Total time,PSNR" > $(prefix)KAPA.csv # CSV header
	@total_time_sum=0; total_time_sq_sum=0; count=0; \
	for i in `seq 1 12`; do \
		output=$$(./sobel_KAPA); \
		time=$$(echo "$$output" | grep "Total time" | awk '{print $$4}'); \
		psnr=$$(echo "$$output" | grep "PSNR" | awk '{print $$9}'); \
		if [ -n "$$time" ]; then \
			echo "$$time,$$psnr" >> $(prefix)KAPA.csv; \
		else \
			echo "Error: could not extract time for execution $$i" >> $(prefix)KAPA.csv; \
		fi; \
	done; \
	python3 calculate_stats.py $(prefix)KAPA
	@echo "Run completed. Results saved in $(prefix)KAPA.csv"
# make stats prefix=test

stats:
	@if [ -z "$(prefix)" ]; then \
		echo "Error: Please provide a 'prefix' argument for the CSV files."; \
		exit 1; \
	fi
	python3 calculate_stats.py $(prefix)

# Makefile rule with parameter for texture17695.bin and CSV output

# Makefile rule with hardcoded input file and CSV output

experiment_12:
	./seq_main -o -b -n 4 -i Image_data/texture17695.bin | grep "Computation timing" | awk '{print $$4}' > texture17695.csv
	./seq_main -o -b -n 4 -i Image_data/texture17695.bin | grep "Computation timing" | awk '{print $$4}' >> texture17695.csv
	./seq_main -o -b -n 4 -i Image_data/texture17695.bin | grep "Computation timing" | awk '{print $$4}' >> texture17695.csv
	./seq_main -o -b -n 4 -i Image_data/texture17695.bin | grep "Computation timing" | awk '{print $$4}' >> texture17695.csv
	./seq_main -o -b -n 4 -i Image_data/texture17695.bin | grep "Computation timing" | awk '{print $$4}' >> texture17695.csv
	./seq_main -o -b -n 4 -i Image_data/texture17695.bin | grep "Computation timing" | awk '{print $$4}' >> texture17695.csv
	./seq_main -o -b -n 4 -i Image_data/texture17695.bin | grep "Computation timing" | awk '{print $$4}' >> texture17695.csv
	./seq_main -o -b -n 4 -i Image_data/texture17695.bin | grep "Computation timing" | awk '{print $$4}' >> texture17695.csv
	./seq_main -o -b -n 4 -i Image_data/texture17695.bin | grep "Computation timing" | awk '{print $$4}' >> texture17695.csv
	./seq_main -o -b -n 4 -i Image_data/texture17695.bin | grep "Computation timing" | awk '{print $$4}' >> texture17695.csv
	./seq_main -o -b -n 4 -i Image_data/texture17695.bin | grep "Computation timing" | awk '{print $$4}' >> texture17695.csv
	./seq_main -o -b -n 4 -i Image_data/texture17695.bin | grep "Computation timing" | awk '{print $$4}' >> texture17695.csv