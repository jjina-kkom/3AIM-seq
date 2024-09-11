#!/bin/bash

############################################################################################################################
##### Input Validation
############################################################################################################################

help(){
        echo "< Usage >"
	echo "          : $0 -f [fastq] -r [read (1/2)] -p [prefix]"
        echo ""
        echo "Arguments:"
        echo "  < Required >"
        echo "     -f or --fastq     Input fastq file name & path"
        echo "     -r or --read      Input fastq read number    (Categorical value: Read 1 -> 1 / Read 2 -> 2)"
	echo "     -p or --prefix    Set prefix"
	echo "  < Optional >" 
        echo "     -o or --output    Define output path         (Default: Working directory)" 
	echo "     -U or --UMI       UMI start position in R2   (Default: 6)"
	echo "     -L or --len       UMI length in R2           (Default: 10)"
        echo "     -h or --help      Print Help (this message) and exit"
}

ARGS=$(getopt -o f:r:o:p:U:L:h --long fastq:,read:,output:,prefix:,UMI:,len:,help -- "$@")

if [[ $? -gt 0 ]]; then
        help
        exit 1
fi

eval set -- ${ARGS}
while :
do
        case $1 in
           -f | --fastq)	FASTQ=$2     ; shift 2 ;;
           -r | --read)		READ=$2      ; shift 2 ;;
           -o | --output)	OUTPUT=$2    ; shift 2 ;;
           -p | --prefix)	PREFIX=$2    ; shift 2 ;;
           -U | --UMI)		UMI_start=$2 ; shift 2 ;;
           -L | --len)		UMI_len=$2   ; shift 2 ;;
           -h | --help)         help         ; exit 0  ;;
           --) shift; break ;;
           *) >&2 echo Unsupported option: $1
                help
                exit 1 ;;
        esac
done

if [[ -z "${FASTQ}" ]] || [[ -z "${READ}" ]] || [[ -z "${PREFIX}" ]]; then
    echo "Error: -f [fastq], -r [read (1/2)], and -p [prefix] are required options."
    help
    exit 1
fi

if [[ -z "${OUTPUT}" ]]    ; then OUTPUT=${PWD}/${PREFIX}; else OUTPUT=${OUTPUT}/${PREFIX}; fi

if [[ -z "${UMI_start}" ]] ; then UMI_start=6                                             ; fi

if [[ -z "${UMI_len}" ]]   ; then UMI_len=10                                              ; fi

############################################################################################################################
##### START 
##############################################################

echo "##### [ Preprocessing from FASTQ (${PREFIX}) ] #####"
echo "=================================================================================================================================================================="
echo "Command :"
if [[ "${READ}" == 2 ]]; then
	echo "Preprocessing_FASTQ.sh --fastq ${FASTQ} --read ${READ} --output ${OUTPUT} --prefix ${PREFIX} --UMI ${UMI_start} --len ${UMI_len}"
else
        echo "Preprocessing_FASTQ.sh --fastq ${FASTQ} --read ${READ} --output ${OUTPUT} --prefix ${PREFIX}"
fi
echo "=================================================================================================================================================================="
echo ""
mkdir -p "${OUTPUT}"

############################################################
start_time=`date '+%Y/%m/%d %H:%M:%S'`
echo "START:	${start_time}"; echo ""

##############################################################
##### 1. FASTQ Seperation 
##############################################################
zcat ${FASTQ} | awk -v O="${OUTPUT}" -v P="${PREFIX}" -v R="${READ}" '
{
	if (NR % 4 == 1) {
		print $0 >> O"/"P"_Header_ID.R"R".txt";
	} else if (NR % 4 == 2) {
		print $0 >> O"/"P"_Sequence.R"R".txt";
	} else if (NR % 4 == 0) { 
		print $0 >> O"/"P"_BaseQ.R"R".txt";
	}
}
'

echo "##### 1. FASTQ Seperation:		`date '+%Y/%m/%d %H:%M:%S'`"

############################################################
##### 2. Split Sequence 
############################################################
awk -v O="${OUTPUT}" -v P="${PREFIX}" -v R="${READ}" '
{
	seq_len = length($0);
	split_seq = "";
	for (i = 1; i <= seq_len; i++) {
		split_seq = split_seq (i == 1 ? "" : "\t") substr($0, i, 1);
	}
	print split_seq >> O"/"P"_Sequence_split.R"R".txt";
}
' ${OUTPUT}/${PREFIX}_Sequence.R${READ}.txt

echo "##### 2. Split Sequence:		`date '+%Y/%m/%d %H:%M:%S'`"

############################################################
##### 3. Convert & Split BaseQ / Calculate stat
############################################################
awk -v O="${OUTPUT}" -v P="${PREFIX}" -v R="${READ}" '
BEGIN {OFS="\t"; 
	for (i = 0; i < 256; i++) {
		ascii[sprintf("%c", i)] = i;
	}
}

{
	sum = 0; 
	sumsq = 0; 
	seq_len = length($0); 
	line = "";

	for (i = 1; i <= seq_len; i++) {
		BQ = ascii[substr($0, i, 1)] - 33; 
		sum += BQ;
		sumsq += BQ * BQ;
		line = line (i == 1 ? "" : "\t") BQ; 
	}

	mean = sum / seq_len;
	variance = (sumsq / seq_len) - (mean * mean);
	stdev = sqrt(variance);

	print line >>  O"/"P"_BaseQ_convert.R"R".txt";
	print mean, stdev >>  O"/"P"_BaseQ_stat.R"R".txt";
}' ${OUTPUT}/${PREFIX}_BaseQ.R${READ}.txt

echo "##### 3. Convert & Split BaseQ:		`date '+%Y/%m/%d %H:%M:%S'`"

############################################################
##### 4. UMI Extraction (Only R2)
############################################################
if [ "${READ}" == "2" ]; then
        awk -v O="${OUTPUT}" -v P="${PREFIX}" -v R="${READ}" -v U="${UMI_start}" -v L="${UMI_len}" '
        {
                UMI = substr($0, U, L);
                print UMI >>  O"/"P"_UMI.txt";
        }' ${OUTPUT}/${PREFIX}_Sequence.R${READ}.txt
	
	cat ${OUTPUT}/${PREFIX}_UMI.txt | sort | uniq -c | awk '{print $2"\t"$1}' | sort -k2 -nr > ${OUTPUT}/${PREFIX}_UMI_summary.txt
	gzip ${OUTPUT}/${PREFIX}_UMI.txt

	echo "##### 4. UMI Extraction (Only R2):	`date '+%Y/%m/%d %H:%M:%S'`"
fi

############################################################
##### 5. Compressing
############################################################

gzip ${OUTPUT}/${PREFIX}_Header_ID.R${READ}.txt
gzip ${OUTPUT}/${PREFIX}_Sequence.R${READ}.txt
gzip ${OUTPUT}/${PREFIX}_Sequence_split.R${READ}.txt
gzip ${OUTPUT}/${PREFIX}_BaseQ.R${READ}.txt
gzip ${OUTPUT}/${PREFIX}_BaseQ_convert.R${READ}.txt

echo "##### 5. Compressing:			`date '+%Y/%m/%d %H:%M:%S'`"

############################################################
end_time=`date '+%Y/%m/%d %H:%M:%S'`
echo ""; echo "FINISH:	${end_time}"
############################################################
