egrep -v '(^C|^c|^!)' *.F90 |\
grep subroutine |\
cut -f1 -d"(" |\
sed -e '/[eE]nd [ ]*subroutine/d' -e 's/:/	/g'  |\
awk '{
   print $3"	"$1"	/" $2" "$3}' > tagsA

egrep -v '(^C|^c|^!)' *.F90 |\
grep module |\
cut -f1 -d"(" |\
sed -e '/[eE]nd [ ]*module/d' -e 's/:/	/g' -e "/procedure/d" |\
awk '{
   print $3"	"$1"	/" $2" "$3}'  > tagsB

export LC_COLLATE=C
cat tagsA tagsB | sort -u --ignore-case > tags
rm tagsA tagsB
