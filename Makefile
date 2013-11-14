include ../c_flags

###########################################################################
#
#   Note on special forms
#      $^ are prereq    $@ is target    $* is stem
#
#   Examples of depending on libaries
#        LDLIBS =  -lthing -lregression 
#        libs_used =  ../lib/libthing.a ../lib/libregression.a 
# 
###########################################################################


PROJECT_NAME = text

# OPT = -O3 -std=c++0x -DNDEBUG

OPT = -O3  -std=c++0x

USES = utils eigen

EXTERNAL_USES = boost_regex gomp

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

level_1 = k_means.o token_manager.o confusion_matrix.o porter.o vocabulary.o eigenword_dictionary.o regex.o
level_2 = helpers.o
level_3 = classifier.o regressor.o
level_4 = bigram.o
level_5 = 

regex_test: regex_test.o
	$(GCC) $^ $(LDLIBS) -o $@
	./regex_test

porter: porter.o
	$(GCC) $^ $(LDLIBS) -o  $@

regressor: regressor.o vocabulary.o regex.o eigenword_dictionary.o helpers.o
	$(GCC) $^ $(LDLIBS) -o  $@

seq_regression: seq_regression.o helpers.o
	$(GCC) $^ $(LDLIBS) -o  $@

lsa: lsa.o vocabulary.o regex.o eigenword_dictionary.o
	$(GCC) $^ $(LDLIBS) -o  $@

bigram: bigram.o k_means.o token_manager.o classifier.o confusion_matrix.o
	$(GCC) $^ $(LDLIBS) -o  $@

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

##################

twain_path = text_src/twain/

##################

tag_path   = text_src/stanford-postagger-2013-04-04/
tag_vers   = -3.1.5

cls_path   = -classpath $(tag_path)stanford-postagger$(tag_vers).jar 
tagger     = edu.stanford.nlp.tagger.maxent.MaxentTagger 
tag_model  = -model $(tag_path)models/wsj-0-18-bidirectional-nodistsim.

# java -classpath stanford-postagger.jar edu.stanford.nlp.process.Morphology
#                -stem  hellos.txt > hello.txt

# script converts to lower case, deletes blank tokens (in call to sed, $$ converts in Make to $)
get_twain: 
	scp sob:/data/gutenberg/twain/*.txt $(twain_path)

# 22 May 2013   6 Twain books: Tagged 982153 words at 3121.78 words per second.
tagged/twain.tagged:
	cat $(twain_path)*.txt | tr '[:upper:]' '[:lower:]' >> tmp.txt
	java -mx2g $(cls_path) $(tagger) $(tag_model) -nthreads 4 -textFile tmp.txt -outputFormat tsv  | sed '/^$$/d' >> tagged/twain.tagged
	rm -f tmp.txt

tagged/ptb45.tagged:
	scp sob:/data/pos_eval_corpora/ptb45/ptb45.tagged  tagged

tagged/ptb17.tagged:
	scp sob:/data/pos_eval_corpora/ptb17/ptb17.tagged  tagged

tagged/train.tagged:
	head -n 1100000 tagged/ptb45.tagged > $@

tagged/test.tagged:
	tail -n 500000 tagged/ptb45.tagged > $@


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  temp directory for all temp text files; emptied by make clean
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

temppath = text_src/temp/

CLEAN = $(temppath)*

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  federalist papers   0 Hamilton, 1 Madison, 11 Ham and Mad, 7 Ham or Mad, 10 Jay
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

fedpath = /Users/bob/data/federalist/

$(temppath)fedregr.txt: $(fedpath)federalist_regr.txt
	tr '[:upper:]' '[:lower:]' < $^ > $@

$(temppath)fedx.txt: $(temppath)fedregr.txt                        # removes leading author
	cut --delimiter=' ' --fields=1 --complement $^ > $@

federalist: regressor lsa $(temppath)fedregr.txt $(temppath)fedx.txt
	./regressor --vocab_file=$(temppath)fedx.txt --regr_file=$(temppath)fedregr.txt  --n_projections 100 --power_iter 1  --bidirectional  
	./lsa --vocab_file=$(temppath)fedx.txt --n_projections 100 --power_iter 1


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  real estate text descriptions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

city = ChicagoOld3

#                  removes lines with no text (need $$ to escape $ in make)
$(temppath)$(city).txt: text_src/real_estate/Set10Tokenized/$(city)Tokenized
	grep -v '^[0-9\,\.[:blank:]]\+$$' $^ > $@

#                  defines random projection
seed  = 2763
#                  number of projections for W, and nProj (each side) for bigram	
nProj = 1500

#                  allowed to have different files for vocabulary and for regression
vFile   = $(temppath)$(city).txt 
rFile   = $(vFile)
outPath = $(temppath)$(city)/

$(outPath)$(nProj).txt: regressor $(vFile) $(rfile) 
	./regressor --vocab_file=$(vFile) --regr_file=$(rFile) --output_path=$(outPath)  -s $(seed) --n_projections $(nProj) --power_iter 1  --bidirectional  
	date >> $@

dore:  $(outPath)$(nProj).txt
	echo 'Running regressor'


# aic sequential regressions

nDocs = 7383

#                          y includes m counts
$(outPath)y.txt: $(outPath)$(nProj).txt
	cut -f 1-2 $(outPath)parsed.txt > $@

$(outPath)aic_%.txt: seq_regression $(outPath)y.txt 
	cut -f 1-$* $(outPath)bigram_$(nProj).txt | ./seq_regression -n $(nDocs) -Y $(outPath)y.txt -X $(outPath)LSA_$(nProj).txt -x $(nProj) -i $* -o $@

doaic: $(outPath)aic_10.txt $(outPath)aic_20.txt $(outPath)aic_30.txt $(outPath)aic_40.txt $(outPath)aic_50.txt $(outPath)aic_75.txt $(outPath)aic_100.txt \
       $(outPath)aic_200.txt $(outPath)aic_400.txt $(outPath)aic_600.txt $(outPath)aic_800.txt $(outPath)aic_900.txt $(outPath)aic_950.txt \
       $(outPath)aic_1000.txt $(outPath)aic_1050.txt $(outPath)aic_1100.txt $(outPath)aic_1200.txt \
       $(outPath)aic_1300.txt $(outPath)aic_1400.txt $(outPath)aic_1500.txt
	echo $(outPath)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  simulated topic model text descriptions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

simpath = text_src/sim/

nSimProj = 200

simFile = $(simpath)sim.txt

$(temppath)sim_regr.txt: regressor $(simFile)
	./regressor --vocab_file=$(simFile) --regr_file=$(simFile) --output_file=$@  -s $(seed) --n_projections $(nSimProj) --power_iter 1  --bidirectional  

dosim:  $(temppath)sim_regr.txt


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  wine descriptions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

winepath = text_src/wine/

nWineProj = 500

vWineFile = $(winepath)RatingsAndNotes.tokens

$(temppath)wine_regr.txt: regressor $(vWineFile)
	./regressor --vocab_file=$(vWineFile) --regr_file=$(vWineFile) --output_file=$@  -s $(seed) --n_projections $(nWineProj) --power_iter 1  --bidirectional  

dowine:  $(temppath)wine_regr.txt


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  google eigenwords
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Error here: cut separates terms with embedded comma in string
epath = text_src/eigenwords/

$(temppath)google.txt: $(epath)pretty_2_grams_PC_100k_300.csv
	sed 's/, / /g' $^ | cut --delimiter=' ' --fields=1,6-26,306-326 > $@

$(epath)pretty_2_grams_PC_100k_300.csv:
	scp sob:/data/pretty/pretty_2_grams_PC_100k_300.csv $(epath)

# unigram vocabulary
vpath = text_src/google/

$(vpath)vocab: $(vpath)vocab.gz
	gunzip $^

$(vpath)vocab.gz:
	scp sob:/data/google_data/1gms/vocab.gz $(vpath)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  classifier application for POS
#   options for folding in other tags, normalizing the bigram rows, weighed avg in clustering, cluster max iterations, tag printing
#  --scale_data is on/off
base_options = --scale_data --bidirectional  --skip 0 --threshold 0.0004 --weight_centroid --iterations 20 --print 0

bigram_test: bigram tagged/train.tagged tagged/test.tagged
	cat tagged/train.tagged | \
	./bigram --projections 100 --clusters 200 --validation tagged/test.tagged  $(base_options)  \
	> results/test/p100_c200

# for debugging, shrink the training/validation data to 100,000 cases and read from std::cin in bigram.cc
bigram_debug: bigram tagged/train.tagged tagged/test.tagged
	./bigram  --bidirectional --projections 100 --clusters 200 --weight_centroid --validation tagged/test.tagged


# test accuracy   bidirectional projections distance clusters
#     0.63            yes           100         2       200    diagonal apparent in confusion; singleton clusters
#     0.27            yes           100        cos      200    clusters appear to be a hodgepodge

# for validation...
#	./bigram                 --projections 200 --distance 2 --clusters 15 --validation tagged/validation.tagged   $(base_options)  \
# for bidirectional...
#	./bigram --bidirectional --projections 200 --distance 2 --clusters 1000 --validation tagged/validation.tagged $(base_options)  \



# ----------------------------------------------------------------------------------------
#  parallel make with fixed number of projections, varying num clusters, both cosine/L2
#  match variables 'task', 'skip', and 'proj' in make command
#         results/$task/skip_$skip/$proj
#
#   make -f -j 4 results/twain/skip_0/050     
#   make -f -j 4 results/ptb45/skip_5/125
#
#  these choices must match the make command
task  = ptb45
skip  =  40
proj  = 125

#  ---  automagic section --- 
tags  =  tagged/$(task).tagged
path  =  results/$(task)/skip_$(skip)/

prefx = $(path)p$(proj)

$(path)/.directory_built: 
	echo Building directory $(path)
	mkdir $(path)
	touch $@

$(path)p$(proj)_dc_c%: bigram $(tags) $(path)/.directory_built
	./bigram  --skip $(skip) --projections $(proj) --distance c --clusters $* $(base_options)  <  $(tags)  >> $@

$(path)p$(proj)_d2_c%: bigram $(tags) $(path)/.directory_built
	./bigram  --skip $(skip) --projections $(proj) --distance 2 --clusters $* $(base_options)  <  $(tags)  >> $@

$(path)$(proj): $(path)/.directory_built \
		$(prefx)_dc_c0050  $(prefx)_dc_c0125  $(prefx)_dc_c0250  $(prefx)_dc_c0500 $(prefx)_dc_c0750 $(prefx)_dc_c1000 \
	        $(prefx)_d2_c0050  $(prefx)_d2_c0125  $(prefx)_d2_c0250  $(prefx)_d2_c0500 $(prefx)_d2_c0750 $(prefx)_d2_c1000
	rm -f $@
	tail -n 1 $(prefx)_dc_c* $(prefx)_d2_c* > $@
	cat $@

###########################################################################
include ../rules_for_makefiles
