include ../c_flags

###########################################################################
#
#   Note on special forms
#      $^ are prereq   $< is first prerequisite  $(word 2,$^) gives the second
#      $@ is target    $* is stem   % -> $*
#
###########################################################################


PROJECT_NAME = text

# OPT = -O3 -std=c++0x -DNDEBUG

OPT = -O3 -fopenmp -std=c++0x

USES = eigen utils

# mpi
EXTERNAL_USES = boost_system boost_thread boost_regex gomp

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

level_1 = k_means.o token_manager.o confusion_matrix.o porter.o vocabulary.o eigenword_dictionary.o 
level_2 = helpers.o
level_3 = classifier.o regressor.o bigram.o unified_regressor.o  lsa_regr.o  bigram_regr.o
level_4 = cluster.o
level_5 = 

regex_test: regex_test.o
	$(GCC) $^ $(LDLIBS) -o $@
	./regex_test

porter: porter.o
	$(GCC) $^ $(LDLIBS) -o  $@

bigram: bigram.o vocabulary.o eigenword_dictionary.o helpers.o
	$(GCC) $^ $(LDLIBS) -o  $@

regressor: regressor.o vocabulary.o regex.o eigenword_dictionary.o helpers.o
	$(GCC) $^ $(LDLIBS) -o  $@

unified_regressor: unified_regressor.o vocabulary.o eigenword_dictionary.o helpers.o
	$(GCC) $^ $(LDLIBS) -o  $@

seq_regression: seq_regression.o helpers.o
	$(GCC) $^ $(LDLIBS) -o  $@

lsa_regr: lsa_regr.o vocabulary.o helpers.o
	$(GCC) $^ $(LDLIBS) -o  $@

bigram_regr: bigram_regr.o vocabulary.o helpers.o
	$(GCC) $^ $(LDLIBS) -o  $@

cluster: cluster.o k_means.o token_manager.o classifier.o confusion_matrix.o
	$(GCC) $^ $(LDLIBS) -o  $@

anes_reply_encoder: anes_reply_encoder.o vocabulary.o eigenword_dictionary.o helpers.o
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
#  movie ratings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

moviePath = text_src/temp/movie_ratings/
movieProj = 500
movieSeed = 28731

$(moviePath)merged.txt: ~/C/tools/merge_movies

$(moviePath)L$(movieProj).txt: lsa_regr $(moviePath)merged_rating_text.txt 
	./lsa_regr --file=$(word 2,$^) --output_path=$(moviePath) -s $(movieSeed) --n_projections $(movieProj) --power_iter 4  --adjustment 'b' --min_frequency 3
	date > $@

domovies: $(moviePath)L$(movieProj).txt
	echo "Done movie projections"


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  real estate text descriptions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

city = ChicagoOld3

#			removes lines with no text (need $$ to escape $ in make)
$(temppath)$(city).txt: text_src/real_estate/Set10Tokenized/$(city)Tokenized
	grep -v '^[0-9\,\.[:blank:]]\+$$' $^ > $@
#			skip for bigram
reSkip = 0
#			defines random projection
reSeed = 2763
#			number of projections for W, and nProj (each side) for bigram	
nProj  = 500
#			allow different files for vocabulary and for regression
vFile   = $(temppath)$(city).txt 
rFile   = $(vFile)
outREPath = $(temppath)$(city)/


engine = regressor

$(outREPath)$(nProj).txt: $(engine) $(vFile) $(rfile) 
	./$(engine) --vocab_file=$(vFile) --regr_file=$(rFile) --output_path=$(outREPath) -s $(reSeed) -k $(reSkip) --n_projections $(nProj) --power_iter 1
	date >> $@

dore:  $(outREPath)$(nProj).txt
	echo 'Running text analysis with engine $(engine)' 


# --- lsa regressions   (dummy targets used to date stamp for linear/quadratic choice in dolsa)
#	adjustments for doc/word matrix are  w(raw, without), t(tfidf), r(row), c(col), n(recip), and b(cca or both)

lsaProj = 1500

#	linear
$(outREPath)L$(lsaProj).txt: lsa_regr $(rfile) 
	./lsa_regr --file=$(rFile) --output_path=$(outREPath) -s $(reSeed) --n_projections $(lsaProj) --power_iter 4  --adjustment 'b' --min_frequency 3
	date >> $@

#	quadratic
$(outREPath)Q$(lsaProj).txt: lsa_regr $(rfile) 
	./lsa_regr --file=$(rFile) --output_path=$(outREPath) -s $(reSeed) --n_projections $(lsaProj) --power_iter 1 --adjustment 'b' --min_frequency 3 --quadratic 
	date >> $@

dolsa:  $(outREPath)L$(lsaProj).txt $(temppath)$(city).txt
	echo '---- LSA analysis ----' 


# --- bigram regressions

bigProj = 1500

$(outREPath)$(bigProj).txt: bigram_regr $(rfile) 
	./bigram_regr --file=$(rFile) --output_path=$(outREPath) -s $(reSeed) --n_projections $(bigProj) --power_iter 4  --adjustment 'b' --min_frequency 3
	date >> $@

dobig:  $(outREPath)$(bigProj).txt $(temppath)$(city).txt
	echo '---- Bigram analysis ----' 


# --- bigram decompositions

bigramPath = $(temppath)bigram/

$(bigramPath).dir_built: 
	mkdir $(bigramPath)
	touch $@

$(bigramPath)date.txt: bigram $(rFile) $(bigramPath).dir_built 
	./bigram --text_file=$(rFile) --output_path=$(bigramPath)  -s $(seed) --n_projections $(nProj) --power_iter 1
	date >> $@

dobigram: $(bigramPath)date.txt


# --- aic sequential regressions

nDocs = 7384

cvseed = 53853

cvPath = $(outPath)cv_$(cvseed)/

prfx = _pre_big_

$(cvPath).directory_built: 
	echo "Building directory for holding cv details."
	mkdir $(cvPath)
	touch $@

#                          y includes m counts
$(outPath)y.txt: $(outPath)$(nProj).txt 
	cut -f 1-2 $(outPath)parsed.txt > $@

#                          first 1500 columns of bigram (left side)
$(outPath)big_$(nProj).txt: $(outPath)bigram_$(nProj).txt 
	cut -f 1-1500 $^ > $@

$(cvPath)aic_pre_big_%.txt: seq_regression $(outPath)y.txt  $(cvPath).directory_built  
	cut -f 1-$* $(outPath)bigram_$(nProj).txt | \
          ./seq_regression -X $(outPath)LSA_$(nProj).txt  -n $(nDocs) -s $(cvseed) -v 10 -Y $(outPath)y.txt -x $(nProj) -i $* -o $@

$(cvPath)aic_pre_lsa_%.txt: seq_regression $(outPath)y.txt  $(cvPath).directory_built
	cut -f 1-$* $(outPath)LSA_$(nProj).txt| \
          ./seq_regression -X $(outPath)big_$(nProj).txt -n $(nDocs) -s $(cvseed) -v 10 -Y $(outPath)y.txt -x $(nProj) -i $* -o $@

doaic: $(cvPath)aic$(prfx)10.txt   $(cvPath)aic$(prfx)20.txt   $(cvPath)aic$(prfx)30.txt   $(cvPath)aic$(prfx)40.txt   $(cvPath)aic$(prfx)50.txt\
       $(cvPath)aic$(prfx)75.txt  $(cvPath)aic$(prfx)100.txt  $(cvPath)aic$(prfx)200.txt  $(cvPath)aic$(prfx)400.txt  $(cvPath)aic$(prfx)600.txt\
      $(cvPath)aic$(prfx)800.txt  $(cvPath)aic$(prfx)900.txt  $(cvPath)aic$(prfx)950.txt $(cvPath)aic$(prfx)1000.txt $(cvPath)aic$(prfx)1050.txt\
     $(cvPath)aic$(prfx)1100.txt $(cvPath)aic$(prfx)1200.txt $(cvPath)aic$(prfx)1300.txt $(cvPath)aic$(prfx)1400.txt $(cvPath)aic$(prfx)1500.txt
	echo $(cvPath)

#                   sets seed=zero to set up testing

$(outPath)big_750.txt:  $(outPath)big_$(nProj).txt
	cut -f 1-750 $^ > $@

testaic: seq_regression $(outPath)y.txt $(outPath)big_750.txt
	cut -f 1-10 $(outPath)LSA_$(nProj).txt | \
        ./seq_regression -X $(outPath)big_750.txt -x 750  -n $(nDocs) -s 0 -v 10 -Y $(outPath)y.txt -i 10 -o $@



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

winePath = text_src/wine/

wineSeed = 36293

nWineProj = 1000

vWineFile = $(winePath)RatingsAndNotes.tokens
outWinePath = $(temppath)wine/


$(outWinePath).directory_built: 
	mkdir $(outWinePath)
	touch $@

$(temppath)wine_regr.txt: regressor $(vWineFile) $(outWinePath).directory_built
	./regressor --vocab_file=$(vWineFile) --regr_file=$(vWineFile) --output_path=$(outWinePath)  -s $(seed) --n_projections $(nWineProj) \
            --power_iter 1  --bidirectional  

$(outWinePath)L$(nWineProj).txt: lsa_regr $(vWineFile) 
	./lsa_regr --file=$(vWineFile) --output_path=$(outWinePath) -s $(wineSeed) --n_projections $(nWineProj) --power_iter 1 --adjustment 'n'
	date >> $@

dowine:  $(outWinePath)L$(nWineProj).txt


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  google eigenwords
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Cut used to pick PCs from left and right context fields
#    See column header to interpret cut columns
#    Error here: cut separates terms with embedded comma in string
# 

epath = text_src/eigenwords/

# trigram eigenwords

$(epath)google_tri.txt: $(epath)pretty_m_3_grams_PHC_50k.csv
	sed 's/, / /g' $^ | cut --delimiter=' ' --fields=1,6-25 > $@

$(epath)pretty_m_3_grams_PHC_50k.csv:
	scp sob:/data/pretty/pretty_m_3_grams_PHC_50k.csv $(epath)

# bigrams

$(temppath)google_bi.txt: $(epath)pretty_2_grams_PC_100k_300.csv
	sed 's/, / /g' $^ | cut --delimiter=' ' --fields=1,6-26,306-326 > $@

$(epath)pretty_2_grams_PC_100k_300.csv:
	scp sob:/data/pretty/pretty_2_grams_PC_100k_300.csv $(epath)

# unigram vocabulary

$(epath)vocab: $(epath)vocab.gz
	gunzip $^

$(epath)vocab.gz:
	scp sob:/data/google_data/1gms/vocab.gz $(epath)

# ANES text

apath = text_src/anes/

$(apath)brown.csv: 
	scp anes.ldc.upenn.edu:/data/anes_revised/coding_dfs/OFCREC.KNOW_BROWN.csv $@

$(apath)cheney.csv: 
	scp anes.ldc.upenn.edu:/data/anes_revised/coding_dfs/OFCREC.KNOW_CHENEY.csv $@

$(apath)pelosi.csv: 
	scp anes.ldc.upenn.edu:/data/anes_revised/coding_dfs/OFCREC.KNOW_PELOSI.csv $@

$(apath)roberts.csv: 
	scp anes.ldc.upenn.edu:/data/anes_revised/coding_dfs/OFCREC.KNOW_ROBERTS.csv $@

# valid names are brown, cheney, pelosi, roberts
name = pelosi

$(apath)$(name).txt: $(apath)$(name).csv $(apath)sed.script
	cut -d ',' -f 2 $< | tail -n +2 | sed -f $(apath)sed.script >> $@

doanes: anes_reply_encoder $(epath)google_tri.txt $(apath)$(name).txt
	./anes_reply_encoder --name $(name)



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
