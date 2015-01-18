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

level_1 = convert.o embed.o
level_2 =
level_3 =

cleanup:
	rm -f prep_events.txt tag_count.txt rectangle_data.txt
	rm -f vocabulary.txt embedded_data.txt

.PHONY: all

all: embedded_data.txt


###########################################################################
#
#	Prepositions
#
#		- first convert into column form
#               - probably want to edit tags.txt to pick subset to keep
#
###########################################################################

nlines = 10000

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

# reverse to overwrite in dictionary with more common term
reversed_eigenwords.en: ~/C/text/eigenwords/eigenwords.300k.200.en.zip
	rm -f $@
	unzip -p $< | tac > $@

embedded_data.txt: rectangle_data.txt vocabulary.txt reversed_eigenwords.en embed
	./embed --eigen_file=reversed_eigenwords.en --eigen_dim 15 --vocab=vocabulary.txt < rectangle_data.txt > $@

###########################################################################

include ~/C/rules_for_makefiles
