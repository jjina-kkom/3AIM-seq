#!/bin/bash

############################################################################################################################
##### Input Validation
############################################################################################################################

help(){
        echo "< Usage >"
	echo "== Without Linker =="
	echo "          : $0 -f [python script] -p [prefix] -P [prefix for filtering info] -l [lenDF file] -r [read list file] --polyA_1 [PolyA length]"
        echo "== With Linker =="
	echo "          : $0 -f [python script] -p [prefix] -P [prefix for filtering info] -l [lenDF file] -r [read list file] --polyA_1 [1st PolyA length] --linker [Linker length] --polyA_2 [2nd PolyA length]"
        echo ""
        echo "Arguments:"
        echo "  < Required >"
	echo "     -f  or --file                 Python script                                                      (3-3_Filtering_Data_byReads.py)"
	echo "     -p  or --samplePrefix         Set prefix for sampleID"
        echo "     -P  or --filterPrefix         Set prefix for filtering step"
        echo "     -l  or --lenDF                DataFrame of measured polyA length (= Output of measurement step)"
        echo "     -r  or --readList             Read list that is survived reads   (= Output of filtering step)"
	echo "            --polyA_1              1st polyA length"
	echo "  < Optional >" 
	echo "     -o  or --output               Directory of output                                                (Default: Working directory)" 
        echo "            --linker               (With Linker) Linker length included in poly A structure"
        echo "            --polyA_2              (With Linker) 2nd polyA length"
        echo "     -h  or --help                 Print Help (this message) and exit"
}

ARGS=$(getopt -o f:p:P:l:r:o:h --long file:,samplePrefix:,filterPrefix:,lenDF:,readList:,output:,help,polyA_1:,linker:,polyA_2: -- "$@")

if [[ $? -gt 0 ]]; then
        help
        exit 1
fi

eval set -- ${ARGS}
while :
do
        case $1 in
	   -f  | --file)              file=$2               ; shift 2 ;;
	   -p  | --samplePrefix)      samplePrefix=$2       ; shift 2 ;;
           -P  | --filterPrefix)      filterPrefix=$2       ; shift 2 ;;
	   -l  | --lenDF)             lenDF=$2              ; shift 2 ;;
	   -r  | --readList)          readList=$2           ; shift 2 ;;
	   -o  | --output)            output=$2             ; shift 2 ;;
                 --polyA_1)           polyA_1=$2            ; shift 2 ;;
	         --linker)            linker=$2             ; shift 2 ;;
                 --polyA_2)           polyA_2=$2            ; shift 2 ;;
           -h  | --help)              help                  ; exit 0  ;;
           --) shift; break ;;
           *) >&2 echo Unsupported option: $1
                help
                exit 1 ;;
        esac
done

if [[ -z "${file}" ]] || [[ -z "${samplePrefix}" ]] || [[ -z "${filterPrefix}" ]] || [[ -z "${lenDF}" ]] || [[ -z "${readList}" ]] || [[ -z "${polyA_1}" ]]; then
    echo "Error: -f, -p, -P, -l, -r and -polyA_1 are required options."
    help
    exit 1
fi

if [[ -z "${output}" ]]   ; then output=${PWD}   ; else output=${output}    ; fi

if [[ -n "${linker}" ]]   ; then linker=$(echo "--Linker_len ${linker}")    ; fi
if [[ -n "${polyA_2}" ]]  ; then polyA_2=$(echo "--polyA_2_len ${polyA_2}") ; fi

############################################################################################################################
##### START 
##############################################################

mkdir -p "${output}"
echo "##### [ Filtered by surviving reads (${samplePrefix}) ] #####"

##@###########################################################
##### 3-3_Filtering_Data_byReads.py
##############################################################
python_command="python3 ${file} --sampleID ${samplePrefix} --filter_prefix ${filterPrefix} --lenDF_file ${lenDF} --read_file ${readList} --output_path ${output} --polyA_1_len ${polyA_1} ${linker} ${polyA_2}"

${python_command}
