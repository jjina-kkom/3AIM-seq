#!/bin/bash

############################################################################################################################
##### Input Validation
############################################################################################################################

help(){
	echo "< Usage >"
        echo "          : $0 -f [python script] -p [prefix] -i [input path (Preprocessing output)] -I [input path (Filtering output)] -c [raed count] -1 [R1 known seq] -2 [R2 known seq]"
        echo ""
        echo "Arguments:"
        echo "  < Required >"
	echo "     -f  or --file                 Python script                                                 (3-2_Selection_SurvivalReads.py)"
	echo "     -p  or --prefix               Set prefix"
	echo "     -i  or --input_1              Directory of input data (= Output of preprocessing step)"
        echo "     -I  or --input_2              Directory of input data (= Output of filtering step)"
        echo "     -c  or --readCount              Total read count of raw data"
	echo "     -1  or --R1_knownSeq          Sequence of R1 known region"
	echo "     -2  or --R2_knownSeq          Sequence of R2 known region"
	echo "  < Optional >" 
	echo "     -o  or --output               Directory of output                                           (Default: Working directory)"
	echo "            --errorRate            NGS error rate                                                (Default: 0.01)"
        echo "            --STdev_filter_whis    A whisker used to filter through standard deviation of baseQ  (Default: 1.5)"
        echo "            --mutCount_filter_whis A whisker used to filter through base mutation count in read  (Default: 1.5)"
        echo "     -h  or --help                 Print Help (this message) and exit"
}

ARGS=$(getopt -o f:p:i:I:c:1:2:o:h --long file:,prefix:,input_1:,input_2:,readCount:,R1_knownSeq:,R2_knownSeq:,output:,help,errorRate:,STdev_filter_whis:,mutCount_filter_whis: -- "$@")

if [[ $? -gt 0 ]]; then
        help
        exit 1
fi

eval set -- ${ARGS}
while :
do
        case $1 in
	   -f  | --file)                 file=$2                 ; shift 2 ;;
	   -p  | --prefix)               prefix=$2               ; shift 2 ;;
           -i  | --input_1)              input_1=$2              ; shift 2 ;;
           -I  | --input_2)              input_2=$2              ; shift 2 ;;
	   -c  | --readCount)            readCount=$2            ; shift 2 ;;
	   -1  | --R1_knownSeq)	         R1_knownSeq=$2          ; shift 2 ;;
	   -2  | --R2_knownSeq)          R2_knownSeq=$2          ; shift 2 ;;
	   -o  | --output)               output=$2               ; shift 2 ;;
	         --errorRate)            errorRate=$2            ; shift 2 ;;
                 --STdev_filter_whis)    STdev_filter_whis=$2    ; shift 2 ;;
	         --mutCount_filter_whis) mutCount_filter_whis=$2 ; shift 2 ;;
           -h  | --help)                 help                    ; exit 0  ;;
           --) shift; break ;;
           *) >&2 echo Unsupported option: $1

                help
                exit 1 ;;
        esac
done

if [[ -z "${file}" ]] || [[ -z "${prefix}" ]] || [[ -z "${input_1}" ]] || [[ -z "${input_2}" ]] || [[ -z "${readCount}" ]] || [[ -z "${R1_knownSeq}" ]] || [[ -z "${R2_knownSeq}" ]]; then
    echo "Error: -f, -p, -i, -I, -c, -1 and -2 are required options."
    help
    exit 1
fi

if [[ -z "${output}" ]]  ; then output=${PWD}  ; else output=${output} ; fi

if [[ -n "${errorRate}" ]]             ; then errorRate=$(echo "--NGS_ErrorRate ${errorRate}")                               ; fi
if [[ -n "${STdev_filter_whis}" ]]     ; then STdev_filter_whis=$(echo "--STdev_filter_whis ${STdev_filter_whis}")           ; fi
if [[ -n "${mutCount_filter_whis}" ]]  ; then mutCount_filter_whis=$(echo "--mutCount_filter_whis ${mutCount_filter_whis}")  ; fi

############################################################################################################################
##### START 
##############################################################

mkdir -p "${output}"
echo "##### [ Selection Survival Reads from Standard Deviaiton of BaseQ & Base Mutaion Count (${prefix}) ] #####"

##############################################################
##### 3-2_Selection_SurvivalReads.py
##############################################################e
python_command="python3 ${file} --sampleID ${prefix} --input_path_1 ${input_1} --input_path_2 ${input_2} --output_path ${output} --read_count ${readCount} --R1_KnownSeq ${R1_knownSeq} --R2_KnownSeq ${R2_knownSeq} ${errorRate} ${STdev_filter_whis} ${mutCount_filter_whis}"

${python_command}
