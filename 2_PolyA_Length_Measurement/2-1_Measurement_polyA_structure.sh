#!/bin/bash

############################################################################################################################
##### Input Validation
############################################################################################################################

help(){
        echo "< Usage >"
	echo "== Without Linker =="
	echo "          : $0 -f [python script] -p [prefix] -i [input path] --polyA_1 [PolyA length] -c [read count] -s [seq size] -1 [R1 known seq] -2 [R2 known seq]"
        echo "== With Linker =="
	echo "          : $0 -f [python script] -p [prefix] -i [input path] --polyA_1 [1st PolyA length] --linker [Linker length] --polyA_2 [2nd PolyA length] -c [read count] -s [seq size] -1 [R1 known seq] -2 [R2 known seq]"
        echo ""
        echo "Arguments:"
        echo "  < Required >"
	echo "     -f  or --file                 Python script                                                      (2-1_Measurement_polyA_structure.py)"
	echo "     -p  or --prefix               Set prefix"
	echo "     -i  or --input                Directory of input data"
	echo "            --polyA_1              1st polyA length"
	echo "     -c  or --readCount            Read count of fastq"
	echo "     -s  or --seqSize              Sequence length of designed amplicon"
	echo "     -1  or --R1_knownSeq          Sequence of R1 known region"
	echo "     -2  or --R2_knownSeq          Sequence of R2 known region"
	echo "  < Optional >" 
	echo "     -o  or --output               Directory of output                                                (Default: Working directory)" 
        echo "            --linker               (With Linker) Linker length included in poly A structure"
        echo "            --polyA_2              (With Linker) 2nd polyA length"
	echo "     -l  or --readLength           Read length                                                        (Default: 301)"
	echo "     -w1 or --windowSize_1         1st window size of sliding window in measured step                 (Default: 5)"
	echo "     -w2 or --windowSize_2         2nd window size of sliding window in measured step                 (Default: 30)"
	echo "     -nc or --neg_cutoff           Minimum number of consecutive negative numbers in measured step    (Default: 3)"
	echo "     -pc or --pos_cutoff           Minimum number of consecutive positive numbers in measured step    (Default: 3)"
	echo "     -n  or --round_n              Decimal placed                                                     (Default: 5)"
	echo "     -S1 or --sliding_chunkSize    Chunk size used during sliding window algorithm execution          (Default: 10000)"
	echo "     -S2 or --measure_chunkSize    Chunk size used during measured step                               (Default: 500)"
	echo "     -m  or --max_workers          Maximum worker of parallel processing process in measured step     (Default: 50)"
        echo "     -h  or --help                 Print Help (this message) and exit"
}

ARGS=$(getopt -o f:p:i:o:c:s:1:2:l:w1:w2:nc:pc:n:S1:S2:m:h --long file:,prefix:,input:,output:,readCount:,seqSize:,R1_knownSeq:,R2_knownSeq:,readLength:,windowSize_1:,windowSize_2:,neg_cutoff:,pos_cutoff:,round_n:,sliding_chunkSize:,measure_chunkSize:,max_workers:,help,polyA_1:,linker:,polyA_2: -- "$@")

if [[ $? -gt 0 ]]; then
        help
        exit 1
fi

eval set -- ${ARGS}
while :
do
        case $1 in
	   -f  | --file)              file=$2               ; shift 2 ;;
	   -p  | --prefix)            prefix=$2             ; shift 2 ;;
           -i  | --input)             input=$2              ; shift 2 ;;
	         --polyA_1)           polyA_1=$2            ; shift 2 ;;
	   -c  | --readCount)         readCount=$2          ; shift 2 ;;
	   -s  | --seqSize)           seqSize=$2            ; shift 2 ;;
	   -1  | --R1_knownSeq)       R1_knownSeq=$2        ; shift 2 ;;
	   -2  | --R2_knownSeq)       R2_knownSeq=$2        ; shift 2 ;;
	   -o  | --output)            output=$2             ; shift 2 ;;
	         --linker)            linker=$2             ; shift 2 ;;
                 --polyA_2)           polyA_2=$2            ; shift 2 ;;
	   -l  | --readLength)        readLength=$2         ; shift 2 ;;
	   -w1 | --windowSize_1)      windowSize_1=$2       ; shift 2 ;;
	   -w2 | --windowSize_2)      windowSize_2=$2       ; shift 2 ;;
	   -nc | --neg_cutoff)        neg_cutoff=$2         ; shift 2 ;;
	   -pc | --pos_cutoff)        pos_cutoff=$2         ; shift 2 ;;
	   -n  | --round_n)           round_n=$2            ; shift 2 ;;
	   -S1 | --sliding_chunkSize) sliding_chunkSize=$2  ; shift 2 ;;
	   -S2 | --measure_chunkSize) measure_chunkSize=$2  ; shift 2 ;;
	   -m  | --max_workers)       max_workers=$2        ; shift 2 ;;
           -h  | --help)              help                  ; exit 0  ;;
           --) shift; break ;;
           *) >&2 echo Unsupported option: $1
                help
                exit 1 ;;
        esac
done

if [[ -z "${file}" ]] || [[ -z "${prefix}" ]] || [[ -z "${input}" ]] || [[ -z "${polyA_1}" ]] || [[ -z "${readCount}" ]] || [[ -z "${seqSize}" ]] || [[ -z "${R1_knownSeq}" ]] || [[ -z "${R2_knownSeq}" ]]; then
    echo "Error: -f, -p, -i, --polyA_1, -c, -s, -1 and -2 are required options."
    help
    exit 1
fi

if [[ -z "${output}" ]]   ; then output=${PWD}   ; else output=${output}  ; fi

if [[ -n "${linker}" ]]           ; then linker=$(echo "--Linker_len ${linker}"); fi
if [[ -n "${polyA_2}" ]]          ; then polyA_2=$(echo "--polyA_2_len ${polyA_2}"); fi
if [[ -n "${readLength}" ]]       ; then readLength=$(echo "--read_length ${readLength}"); fi
if [[ -n "${windowSize_1}" ]]     ; then windowSize_1=$(echo "--window_size_1 ${windowSize_1}"); fi
if [[ -n "${windowSize_2}" ]]     ; then windowSize_2=$(echo "--window_size_2 ${windowSize_2}"); fi
if [[ -n "${neg_cutoff}" ]]       ; then neg_cutoff=$(echo "--neg_cutoff ${neg_cutoff}"); fi
if [[ -n "${pos_cutoff}" ]]       ; then pos_cutoff=$(echo "--pos_cutoff ${pos_cutoff}"); fi
if [[ -n "${round_n}" ]]          ; then round_n=$(echo "--round_n ${round_n}"); fi
if [[ -n "${sliding_chunkSize}" ]]; then sliding_chunkSize=$(echo "--sliding_chunkSize ${sliding_chunkSize}"); fi
if [[ -n "${measure_chunkSize}" ]]; then measure_chunkSize=$(echo "--measure_chunkSize ${measure_chunkSize}"); fi
if [[ -n "${max_workers}" ]]      ; then max_workers=$(echo "--max_workers ${max_workers}"); fi

############################################################################################################################
##### START 
##############################################################

mkdir -p "${output}"
echo "##### [ Measurement PolyA Length (${prefix}) ] #####"

##############################################################
##### 2-1_Measurement_polyA_structure.py 
##############################################################
python_command="python3 ${file} --sampleID ${prefix} --input_path ${input} --output_path ${output} --polyA_1_len ${polyA_1} ${linker} ${polyA_2} --read_count ${readCount} --seq_size ${seqSize} --R1_KnownSeq ${R1_knownSeq} --R2_KnownSeq ${R2_knownSeq} ${readlength} ${windowSize_1} ${windowSize_2} ${neg_cutoff} ${pos_cutoff} ${round_n} ${sliding_chunkSize} ${measure_chunkSize} ${max_workers}"

${python_command}
