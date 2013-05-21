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

USES = utils

level_1 = k_means.o token_manager.o
level_2 = bigram.o
level_3 = 

##################

text_path = text_src/twain/

tag_path   = stanford-postagger-2013-04-04/
tag_vers   = -3.1.5

cls_path   = -classpath $(tag_path)stanford-postagger$(tag_vers).jar 
tagger     = edu.stanford.nlp.tagger.maxent.MaxentTagger 
tag_model  = -model $(tag_path)models/wsj-0-18-bidirectional-nodistsim.tagger

##################

clean_txt:
	rm -f tokens.txt tmp.txt

# script converts to lower case, deletes blank tokens (in call to sed, $$ converts in Make to $)

tokens.txt:
	cat $(text_path)*.txt | tr '[:upper:]' '[:lower:]' >> tmp.txt
	java -mx2g $(cls_path) $(tagger) $(tag_model) -nthreads 4 -textFile tmp.txt -outputFormat tsv  | sed '/^$$/d' >> tokens.txt

# compute svd of random projected bigram matrix
bigram.o: bigram.cc

bigram: bigram.o k_means.o token_manager.o
	$(GCC) $^ $(LDLIBS) -o  $@

bigram.prj: tokens.txt bigram
	./bigram --threshold 0.005 --projections 100 --scaling 0 --clusters 60 --iterations 20 --print 5 < tokens.txt >> bigram.prj



###########################################################################
include ../rules_for_makefiles
