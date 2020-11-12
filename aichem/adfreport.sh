dirname="$(dirname $1)"
jobname="$(dirname $dirname)"

basename="$(basename $1)"
moleculename="${basename%%.*}"

jobnumber="${basename%.*}"
number="${jobnumber#*.}"

frag1name="${jobname}"/frag1."${number}"/frag1."${number}".t21
frag2name="${jobname}"/frag2."${number}"/frag2."${number}".t21

figurepath="$(dirname $5)"
frag1figurepath="${figurepath}"/"${moleculename}"."${number}"_frag1.html
compfigurepath="${figurepath}"/"${moleculename}"."${number}"_complex.html
frag2figurepath="${figurepath}"/"${moleculename}"."${number}"_frag2.html


frag1=$2
if [ "${frag1[0]}" != 'None' ]; then

  temp2=`echo "$2" | tr -d '"'| tr -d "'"`
  adfreport -i "${frag1name}" ${temp2} -v "-grid Fine" -v "-antialias"  -v "-viewplane {1 0 0}" -v "-bgcolor #ffffff" -o "${frag1figurepath}"

fi

frag2=$4
if [ "${frag2[0]}" != 'None' ]; then

  temp4=`echo "$4" | tr -d '"'| tr -d "'"`
  adfreport -i "${frag2name}" ${temp4} -v "-grid Fine" -v "-antialias"  -v "-viewplane {1 0 0}" -v "-bgcolor #ffffff" -o "${frag2figurepath}"

fi


temp3=`echo "$3" | tr -d '"'| tr -d "'"`
adfreport -i $1 ${temp3} -v "-grid Fine" -v "-antialias"  -v "-viewplane {1 0 0}" -v "-bgcolor #ffffff" -o "${compfigurepath}"
