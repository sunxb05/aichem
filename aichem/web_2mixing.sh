#!/bin/sh

dirname="$(dirname $1)"
jobname="$(dirname $dirname)"

basename="$(basename $1)"
moleculename="${basename%%.*}"

jobnumber="${basename%.*}"
number="${jobnumber#*.}"

figurepath="$(dirname $8)"
finalfigurepath="${figurepath}"/"${moleculename}"."${number}".html

frag1figurepath="${moleculename}"."${number}"_frag1
compfigurepath="${moleculename}"."${number}"_complex
frag2figurepath="${moleculename}"."${number}"_frag2

frag1figurepath_full="${figurepath}"/"${moleculename}"."${number}"_frag1.html
compfigurepath_full="${figurepath}"/"${moleculename}"."${number}"_complex.html
frag2figurepath_full="${figurepath}"/"${moleculename}"."${number}"_frag2.html

informpath="${figurepath}"/"${moleculename}"."${number}"_inform.html
adfreport -i $1 'Molecule' 'Symmetry'  "BondingEnergy"  "Geometry%xyz#12.4f##3" -v "-grid Fine" -v "-antialias"  -v "-viewplane {1 0 0}" -v "-bgcolor #ffffff" -o "${informpath}"

declare -a frag1orbitallist
declare -a complexorbtiallist
declare -a frag2orbitallist

frag1orbitallist[0]=$2
frag1orbitallist[1]=$3
complexorbtiallist[0]=$4
complexorbtiallist[1]=$5
frag2orbitallist[0]=$6
frag2orbitallist[1]=$7

for orbital in "${complexorbtiallist[0]}"
  do
    IFS=";" read -r -a arr <<< "${orbital}"
    mix1_comp_num=${#arr[@]}
    if [ ${mix1_comp_num} == 3 ]; then
      # mix1_comp_1=${arr[0]}
      mix1_comp_1=`echo "${arr[0]}" | tr -d "'"`
      mix1_comp_2=`echo "${arr[1]}" | tr -d "'"`
      mix1_comp_3=`echo "${arr[2]}" | tr -d "'"`

      mix1_compf_1=0
      mix1_compf_2=1
      mix1_compf_3=2
    else
      # mix1_comp_1=${arr[0]}
      mix1_comp_1=`echo "${arr[0]}" | tr -d "'"`
      mix1_comp_2=' '
      mix1_comp_3=`echo "${arr[1]}" | tr -d "'"`

      mix1_compf_1=0
      mix1_compf_2=' '
      mix1_compf_3=1
    fi
  done

for orbital in "${complexorbtiallist[1]}"
  do
    IFS=";" read -r -a arr <<< "${orbital}"
    mix2_comp_num=${#arr[@]}
    if [ ${mix2_comp_num} == 3 ]; then
      mix2_comp_1=`echo "${arr[0]}" | tr -d "'"`
      mix2_comp_2=`echo "${arr[1]}" | tr -d "'"`
      mix2_comp_3=`echo "${arr[2]}" | tr -d "'"`

      mix2_compf_1=$((0 + mix1_comp_num))
      mix2_compf_2=$((1 + mix1_comp_num))
      mix2_compf_3=$((2 + mix1_comp_num))
    else
      mix2_comp_1=`echo "${arr[0]}" | tr -d "'"`
      mix2_comp_2=' '
      mix2_comp_3=`echo "${arr[1]}" | tr -d "'"`

      mix2_compf_1=$((0 + mix1_comp_num))
      mix2_compf_2=' '
      mix2_compf_3=$((1 + mix1_comp_num))
    fi
  done

for orbital in "${frag1orbitallist[0]}"
  do
    IFS=";" read -r -a arr <<< "${orbital}"
    mix1_frag1_num=${#arr[@]}
    if [ ${mix1_frag1_num} == 2 ]; then
      mix1_frag1_1=`echo "${arr[0]}" | tr -d "'"`
      mix1_frag1_2=`echo "${arr[1]}" | tr -d "'"`

      mix1_frag1f_1=0
      mix1_frag1f_2=1
    else
      mix1_frag1_1=`echo "${arr[0]}" | tr -d "'"`
      mix1_frag1_2=' '

      mix1_frag1f_1=0
      mix1_frag1f_2=' '
    fi
  done

for orbital in "${frag1orbitallist[1]}"
  do
    IFS=";" read -r -a arr <<< "${orbital}"
    mix2_frag1_num=${#arr[@]}
    if [ ${mix2_frag1_num} == 2 ]; then
      mix2_frag1_1=`echo "${arr[0]}" | tr -d "'"`
      mix2_frag1_2=`echo "${arr[1]}" | tr -d "'"`

      mix2_frag1f_1=$((0 + mix1_frag1_num))
      mix2_frag1f_2=$((1 + mix1_frag1_num))
    else
      mix2_frag1_1=`echo "${arr[0]}" | tr -d "'"`
      mix2_frag1_2=' '

      mix2_frag1f_1=$((0 + mix1_frag1_num))
      mix2_frag1f_2=' '
    fi
  done

for orbital in "${frag2orbitallist[0]}"
  do
    IFS=";" read -r -a arr <<< "${orbital}"
    mix1_frag2_num=${#arr[@]}
    if [ ${mix1_frag2_num} == 2 ]; then
      mix1_frag2_1=`echo "${arr[0]}" | tr -d "'"`
      mix1_frag2_2=`echo "${arr[1]}" | tr -d "'"`

      mix1_frag2f_1=0
      mix1_frag2f_2=1
    else
      mix1_frag2_1=`echo "${arr[0]}" | tr -d "'"`
      mix1_frag2_2=' '

      mix1_frag2f_1=0
      mix1_frag2f_2=' '

    fi
  done

for orbital in "${frag2orbitallist[1]}"
  do
    IFS=";" read -r -a arr <<< "${orbital}"
    mix2_frag1_num=${#arr[@]}
    if [ ${mix2_frag1_num} == 2 ]; then
      mix2_frag2_1=`echo "${arr[0]}" | tr -d "'"`
      mix2_frag2_2=`echo "${arr[1]}" | tr -d "'"`

      mix2_frag2f_1=$((0 + mix1_frag2_num))
      mix2_frag2f_2=$((1 + mix1_frag2_num))
    else
      mix2_frag2_1=`echo "${arr[0]}" | tr -d "'"`
      mix2_frag2_2=' '

      mix2_frag2f_1=$((0 + mix1_frag2_num))
      mix2_frag2f_2=' '

    fi
  done

SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"
htmlindex="$SCRIPTPATH"/index.html
cp $htmlindex  $finalfigurepath



echo '<tr>'                                         >> $finalfigurepath
echo '<th> </th>'                                   >> $finalfigurepath
echo '<td class="top"> </a></td>'                   >> $finalfigurepath
echo '<th>'$mix1_comp_1'</th>'                      >> $finalfigurepath
echo '<td class="top"><img src="'$compfigurepath'.jpgs/'$mix1_compf_1'.jpg" width="200" height="200" onerror="this.style.display='\'none\''"></a></td>'    >> $finalfigurepath
echo '<th> </th>'                                   >> $finalfigurepath
echo '<td class="top"> </a></td>'                   >> $finalfigurepath

echo '<th> </th>'                                   >> $finalfigurepath
echo '<td class="top"> </a></td>'                   >> $finalfigurepath
echo '<th> </th>'                                   >> $finalfigurepath
echo '<td class="top"> </a></td>'                   >> $finalfigurepath
echo '<th> </th>'                                   >> $finalfigurepath
echo '<td class="top"> </a></td>'                   >> $finalfigurepath

echo '<th> </th>'                                   >> $finalfigurepath
echo '<td class="top"> </a></td>'                   >> $finalfigurepath
echo '<th>'$mix2_comp_1'</th>'                      >> $finalfigurepath
echo '<td class="top"><img src="'$compfigurepath'.jpgs/'$mix2_compf_1'.jpg" width="200" height="200" onerror="this.style.display='\'none\''"></a></td>'    >> $finalfigurepath
echo '<th> </th>'                                   >> $finalfigurepath
echo '<td class="top"> </a></td>'                   >> $finalfigurepath
echo '</tr>'                                        >> $finalfigurepath




echo '<tr>'                                         >> $finalfigurepath
echo '<th>'$mix1_frag1_1'</th>'                     >> $finalfigurepath
echo '<td class="top"><img src="'$frag1figurepath'.jpgs/'$mix1_frag1f_1'.jpg" width="200" height="200" onerror="this.style.display='\'none\''"></a></td>'    >> $finalfigurepath
echo '<th>'$mix1_comp_2'</th>'                      >> $finalfigurepath
echo '<td class="top"><img src="'$compfigurepath'.jpgs/'$mix1_compf_2'.jpg" width="200" height="200" onerror="this.style.display='\'none\''"></a></td>'    >> $finalfigurepath
echo '<th>'$mix1_frag2_1'</th>'                     >> $finalfigurepath
echo '<td class="top"><img src="'$frag2figurepath'.jpgs/'$mix1_frag2f_1'.jpg" width="200" height="200" onerror="this.style.display='\'none\''"></a></td>'    >> $finalfigurepath

echo '<th> </th>'                                   >> $finalfigurepath
echo '<td class="top"> </a></td>'                   >> $finalfigurepath
echo '<th> </th>'                                   >> $finalfigurepath
echo '<td class="top"> </a></td>'                   >> $finalfigurepath
echo '<th> </th>'                                   >> $finalfigurepath
echo '<td class="top"> </a></td>'                   >> $finalfigurepath


echo '<th>'$mix2_frag1_1'</th>'                     >> $finalfigurepath
echo '<td class="top"><img src="'$frag1figurepath'.jpgs/'$mix2_frag1f_1'.jpg" width="200" height="200" onerror="this.style.display='\'none\''"></a></td>'    >> $finalfigurepath
echo '<th>'$mix2_comp_2'</th>'                      >> $finalfigurepath
echo '<td class="top"><img src="'$compfigurepath'.jpgs/'$mix2_compf_2'.jpg" width="200" height="200" onerror="this.style.display='\'none\''"></a></td>'    >> $finalfigurepath
echo '<th>'$mix2_frag2_1'</th>'                     >> $finalfigurepath
echo '<td class="top"><img src="'$frag2figurepath'.jpgs/'$mix2_frag2f_1'.jpg" width="200" height="200" onerror="this.style.display='\'none\''"></a></td>'    >> $finalfigurepath
echo '</tr>'                                        >> $finalfigurepath



echo '<tr>'                                         >> $finalfigurepath
echo '<th>'$mix1_frag1_2'</th>'                     >> $finalfigurepath
echo '<td class="top"><img src="'$frag1figurepath'.jpgs/'$mix1_frag1f_2'.jpg" width="200" height="200" onerror="this.style.display='\'none\''"></a></td>'    >> $finalfigurepath
echo '<th>'$mix1_comp_3'</th>'                      >> $finalfigurepath
echo '<td class="top"><img src="'$compfigurepath'.jpgs/'$mix1_compf_3'.jpg" width="200" height="200" onerror="this.style.display='\'none\''"></a></td>'    >> $finalfigurepath
echo '<th>'$mix1_frag2_2'</th>'                     >> $finalfigurepath
echo '<td class="top"><img src="'$frag2figurepath'.jpgs/'$mix1_frag2f_2'.jpg" width="200" height="200" onerror="this.style.display='\'none\''"></a></td>'    >> $finalfigurepath

echo '<th> </th>'                                   >> $finalfigurepath
echo '<td class="top"> </a></td>'                   >> $finalfigurepath
echo '<th> </th>'                                   >> $finalfigurepath
echo '<td class="top"> </a></td>'                   >> $finalfigurepath
echo '<th> </th>'                                   >> $finalfigurepath
echo '<td class="top"> </a></td>'                   >> $finalfigurepath


echo '<th>'$mix2_frag1_2'</th>'                     >> $finalfigurepath
echo '<td class="top"><img src="'$frag1figurepath'.jpgs/'$mix2_frag1f_2'.jpg" width="200" height="200" onerror="this.style.display='\'none\''"></a></td>'    >> $finalfigurepath
echo '<th>'$mix2_comp_3'</th>'                      >> $finalfigurepath
echo '<td class="top"><img src="'$compfigurepath'.jpgs/'$mix2_compf_3'.jpg" width="200" height="200" onerror="this.style.display='\'none\''"></a></td>'    >> $finalfigurepath
echo '<th>'$mix2_frag2_2'</th>'                     >> $finalfigurepath
echo '<td class="top"><img src="'$frag2figurepath'.jpgs/'$mix2_frag2f_2'.jpg" width="200" height="200" onerror="this.style.display='\'none\''"></a></td>'    >> $finalfigurepath
echo '</tr>'                                        >> $finalfigurepath

echo '</table></body></html>'                       >> $finalfigurepath

rm $frag1figurepath_full $compfigurepath_full $frag2figurepath_full

open  $informpath $finalfigurepath

