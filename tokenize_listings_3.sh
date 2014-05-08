#!/bin/sh
#
# Each input  line is asking-price followed by the text of the description.
cat $1 |
sed '
s/,/ , /g
s/[.]/ . /g
s/!/ ! /g
s/[(]/ ( /g
s/[)]/ ) /g
s/;/ ; /g
s/:/ : /g
s/\// \/ /g
s/[*]/ * /g
s/["]/ " /g
s/~~~/ /g
s/  */ /g' |
tr A-Z a-z |
sed '
s/ w \/ / w\/ /g
s/\([0-9]\) [.] \([0-9]\)/\1.\2/g
s/\([0-9]\) , \([0-9]\)/\1,\2/g
s/\([a-z]\)- /\1 - /g
s/ 1 \/ 2 / 1\/2 /g
s/ [.] [.] [.] / ... /g' | 
sed -f tokenize_quotes.sed 

#
