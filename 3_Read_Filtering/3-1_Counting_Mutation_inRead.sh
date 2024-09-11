#!/bin/bash

############################################################################################################################
##### Input Validation
############################################################################################################################

help(){
	echo "< Usage >"
        echo "          : $0 -f [python script] -p [prefix] -i [input path] -1 [R1 known seq] -2 [R2 known seq]"
        echo ""
        echo "Arguments:"
        echo "  < Required >"
	echo "     -f  or --file                 Python script                          (3-1_Counting_Mutation_inRead.py)"
	echo "     -p  or --prefix               Set prefix"
	echo "     -i  or --input                Directory of input data"
	echo "     -1  or --R1_knownSeq          Sequence of R1 known region"
	echo "     -2  or --R2_knownSeq          Sequence of R2 known region"
	echo "  < Optional >" 
	echo "     -o  or --output               Directory of output                    (Default: Working directory)" 
	echo "            --R1_knownSeq_start    Start point of known region in R1      (Default: 1)"
        echo "            --R2_knownSeq_start    Start point of known region in R2      (Default: 1)"
        echo "     -h  or --help                 Print Help (this message) and exit"
}

ARGS=$(getopt -o f:p:i:1:2:o:h --long file:,prefix:,input:,R1_knownSeq:,R2_knownSeq:,output:,help,R1_knownSeq_start:,R2_knownSeq_start: -- "$@")

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
           -i  | --input)                input=$2                ; shift 2 ;;
	   -1  | --R1_knownSeq)	         R1_knownSeq=$2          ; shift 2 ;;
	   -2  | --R2_knownSeq)          R2_knownSeq=$2          ; shift 2 ;;
	   -o  | --output)               output=$2               ; shift 2 ;;
	         --R1_knownSeq_start)    R1_knownSeq_start=$2    ; shift 2 ;;
                 --R2_knownSeq_start)    R2_knownSeq_start=$2    ; shift 2 ;;
           -h  | --help)                 help                    ; exit 0  ;;
           --) shift; break ;;
           *) >&2 echo Unsupported option: $1

                help
                exit 1 ;;
        esac
done

if [[ -z "${file}" ]] || [[ -z "${prefix}" ]] || [[ -z "${input}" ]] || [[ -z "${R1_knownSeq}" ]] || [[ -z "${R2_knownSeq}" ]]; then
    echo "Error: -f, -p, -i, -1 and -2 are required options."
    help
    exit 1
fi

if [[ -z "${output}" ]]   ; then output=${PWD}   ; else output=${output}  ; fi

if [[ -n "${R1_knownSeq_start}" ]]  ; then R1_knownSeq_start=$(echo "--R1_KnownSeq_start ${R1_knownSeq_start}")  ; fi
if [[ -n "${R2_knownSeq_start}" ]]  ; then R2_knownSeq_start=$(echo "--R2_KnownSeq_start ${R2_knownSeq_start}")  ; fi

############################################################################################################################
##### START 
##############################################################

mkdir -p "${output}"
echo "##### [ Counting Mutation in Read (${prefix}) ] #####"

##############################################################
##### 3-1_Counting_Mutation_inRead.py
##############################################################
python_command="python3 ${file} --sampleID ${prefix} --input_path ${input} --output_path ${output} --R1_KnownSeq ${R1_knownSeq} --R2_KnownSeq ${R2_knownSeq} ${R1_knownSeq_start} ${R2_knownSeq_start}"

${python_command}

