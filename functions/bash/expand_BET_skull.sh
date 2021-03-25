#!/usr/bin/env bash



HELP() {
    cat <<HELP
expand_BET_skull

Usage: 
    bash ${0##*/} -i <anat-file> -skull <skull-mask> -iskull <inskull-mask> -os <outskin-mask> [options]
    
Compulsory arguments:
    -i            Anatomic image (usually T1w)
    -skull        Skull mask provided by BET
    -iskull       Inskull mask provided by BET
    -skin         Outskin mask provided by BET
    
Optional arguments:
    -os <suffix>  Suffix for the expanded skull image (default is "_exp-skull")
    -n <n>        Nth percentile to use as threshold for skull intensity (default=50)
    -k            Will keep temporary files.
    -p <p>        Prefix for running FSL functions (can be a path or just a prefix)
    
HELP
}

#/********************/
#/********************/

if [[ $# -eq 0 ]]; then
  HELP >&2
  exit 1
fi

# Defaults

FSLPREFIX=""
NP=50
KTMP="No"
OUT_SUFFIX="_exp-skull"



ARGS="$@"
SYN=">> Incorrect syntax. See -help"

while [ "$#" -gt 0 ]
do
  case "$1" in

    -i | -I )
    if [[ -n "$2" ]]; then IN="$2"; else echo $SYN; exit 1 ; fi
    shift 2
    ;;
    -skull |-SKULL )
    if [[ -n "$2" ]]; then SKULL="$2"; else echo $SYN; exit 1 ; fi
    shift 2
    ;;
    -iskull |-ISKULL )
    if [[ -n "$2" ]]; then INSKULL="$2"; else echo $SYN; exit 1 ; fi
    shift 2
    ;;
    -skin |-SKIN )
    if [[ -n "$2" ]]; then SKIN="$2"; else echo $SYN; exit 1 ; fi
    shift 2
    ;;
    -os |-OS )
    if [[ -n "$2" ]]; then OUT_SUFFIX="$2"; else echo $SYN; exit 1 ; fi
    shift 2
    ;;
    -n |-N )
    if [[ -n "$2" ]]; then NP="$2"; else echo $SYN; exit 1 ; fi
    shift 2
    ;;
    -k |-K )
    KTMP="Yes"
    shift
    ;;
    -p |-P )
    if [[ -n "$2" ]]; then FSLPREFIX="$2"; else echo $SYN; exit 1 ; fi
    shift 2
    ;;
    -help|--help)
    HELP >&2
    exit 1
    shift
    ;;
    *)
    if [[ $arguments -eq 0 ]] && [[ "$1" != "" ]]
    then
    echo "${0##*/} >> Unknown argument '$1'. See -help." >&2;
    exit 1
    fi
    shift 1
    ;;
  esac
done






# **************************
# Functions

extract_base_name() {
IN=$1
INext=${IN##*.}
local INbase=""
if [[ $INext == "gz" ]]; then
  INngz=${IN%.gz}
  INext=".${INngz##*.}.gz"
  INbase=${IN%$INext}
else
  INext=".${INext}"
  INbase=${IN%$INext}  
fi
INbase=${INbase##*/}
echo $INbase
}

extract_path() {
IN=$1
BASE=${IN##*/}  # base
local DIR=${IN%$BASE}  # dirpath
echo $DIR
}

# *************************




# *************************
# Files
INname=`extract_base_name $IN`
INpath=`extract_path $IN`

OUT_PREFIX="${INpath}${INname}${OUT_SUFFIX}"


# **************************






cat <<REPORTPARAMETERS

--------------------------------------------------------------------------------------
 Parameters
--------------------------------------------------------------------------------------

 Input file:         $IN
 Skull mask:         $SKULL
 Inskull mask:       $INSKULL
 Outskin mask:       $SKIN
 Output prefix:      $OUT_PREFIX
 Nth percentile:     $NP
 Keep tmp files:     $KTMP
 FSL prefix:         $FSLPREFIX
======================================================================================

REPORTPARAMETERS



# Temp files
if [[ "$OSTYPE" == "darwin" ]]; then # if macOSX
  declare -a TMP # temporary files
else
  declare -A TMP # temporary files
fi

TMP[INskull]="${INpath}${INname}_tmpskull.nii.gz"
TMP[INskin]="${INpath}${INname}_tmpskin.nii.gz"







# Get mean skull intensity
echo "Computing skull threshold..."
"${FSLPREFIX}fslmaths" $IN -mas $SKULL ${TMP[INskull]}
skull_int=`"${FSLPREFIX}fslstats" ${TMP[INskull]} -P $NP` # get smallest roi of non zero voxels
echo "Skull threshold: ${skull_int}"

echo "Thresholding input image"
"${FSLPREFIX}fslmaths" $IN -mas $SKIN ${TMP[INskin]}
"${FSLPREFIX}fslmaths" ${TMP[INskin]} -uthr $skull_int -bin -sub $INSKULL -add $SKULL -bin $OUT_PREFIX

# Remove all temporary files
if [[ $KTMP == "No" ]]; then
  for key in ${!TMP[@]}; do
    if test -f "${TMP[$key]}"; then
      rm ${TMP[$key]}
    fi
  done
fi


exit 0




