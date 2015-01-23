include ~/C/c_flags

###########################################################################
#
#   Note on special forms
#      $^ are prereq   $< is first prerequisite  $(word 2,$^) gives the second
#      $@ is target    $* is stem   % -> $*
#
###########################################################################


PROJECT_NAME = prep_error

OPT = -O3 -std=c++11

# OPT = - -DNDEBUG O3

level_1 = convert.o embed.o embed_auction.o recode_data.o
level_2 =
level_3 =

cleanup:
	rm -f prep_events.txt tag_count.txt rectangle_data.txt
	rm -f vocabulary.txt embedded_data.txt
	rm -f reversed_eigenwords.en 
	rm -rf auction_data

.PHONY: all

all: auction_data


###########################################################################
#
#	Prepositions
#
#		- first convert into column form
#               - probably want to edit tags.txt to pick subset to keep
#
###########################################################################

nlines = 1000

prep_events.txt: ~/data/joel/7m-4d-Aug30-events.gz
	echo Building base file 'prep_events.txt' from $(nlines) sentences.
	gunzip -c $< | head -n $(nlines) > $@

# edit these by hand to define tags.txt  (NINw will add = pairs)
all_tags.txt : prep_events.txt   # delete everything *except* tags ... uniq -c produces mult space in output... damn
	rm -f tags.txt tag_count.txt
	sed -e 's/.*DY //' -e "s/#[-#A-Za-z0-9_?,=!;:\`\_\.\'$$]* / /g" -e 's/ $$//' $< | tr ' ' '\n' | sort | uniq > tag_count.txt
	cat tag_count.txt
	cut -f2,2 -d ' ' tag_count.txt > $@

convert: convert.o
	$(GCC) $^ $(LDLIBS) -o  $@

rectangle_data.txt: prep_events.txt tags.txt convert
	./convert --tag_file=tags.txt < prep_events.txt > $@
	head $@

vocabulary.txt: rectangle_data.txt Makefile   # wipe out header line at start, blank at end
	tail -n +2 $< | tr '\t' '\n' | tr '=' '\n' | sort | uniq | tail -n +2 > $@

embed: embed.o
	$(GCC) $^ $(LDLIBS) -o  $@

embed_auction: embed_auction.o
	$(GCC) $^ $(LDLIBS) -o  $@

recode_data: recode_data.o
	$(GCC) $^ $(LDLIBS) -o  $@

build_cv: build_cv.o
	$(GCC) $^ $(LDLIBS) -o  $@

# reverse Zipf order so later, more common coord overwrite in dictionary
reversed_eigenwords.en: ~/data/text/eigenwords/eigenwords.300k.200.en.gz
	rm -f $@
	gunzip -c $< | tac > $@

embedded_data.txt: rectangle_data.txt vocabulary.txt reversed_eigenwords.en embed
	./embed --eigen_file=reversed_eigenwords.en --eigen_dim 15 --vocab=vocabulary.txt < rectangle_data.txt > $@


jason_dat: csv_parser data/jason_97.csv
	rm -rf data/jason.tar; mkdir data/jason.tar
	./csv_parser -f data/jason_97.csv -d data/jason.tar
	cd data/jason.tar; chmod +x index.sh; ./index.sh > jason.dat

jason_test: auction.test data/jason.tar/index.sh
	 ./auction.test -f data/jason.tar/jason.dat -o data/jason.tar/ -r 100 -a 2 -p 3 -c 5


auction_data: rectangle_data.txt vocabulary.txt reversed_eigenwords.en embed_auction
	rm -rf $@
	mkdir auction_data
	./embed_auction --eigen_file=reversed_eigenwords.en --eigen_dim 5 --vocab=vocabulary.txt -o $@ < rectangle_data.txt
	chmod +x $@/index.sh

word0  = in
word1  = to
inDir  = auction_data
outDir = $(inDir)/$(word0)_$(word1)

.PHONY: build_auction_data run_auction

build_auction_data: recode_data $(inDir)
	rm -rf $(outDir); mkdir $(outDir)
	cp $(inDir)/index.sh $(outDir)
	./recode_data --input_dir=$(inDir) --output_dir=$(outDir) --word0=$(word0) --word1=$(word1)

run_auction: build_auction_data recode_data $(auction_data)
	./auction



###########################################################################

include ~/C/rules_for_makefiles
