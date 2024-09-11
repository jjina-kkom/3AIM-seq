#!/usr/bin/env python
# coding: utf-8

##########################
##### Import modules #####
##########################
import os
import time
import argparse
import numpy as np
import pandas as pd
from alive_progress import alive_bar
from concurrent.futures import ProcessPoolExecutor, as_completed

####################
##### argument #####
####################
parser = argparse.ArgumentParser(description='Measurement poly(A) structure')

parser.add_argument('--sampleID')
parser.add_argument('--input_path')
parser.add_argument('--output_path', default=os.getcwd())
parser.add_argument('--polyA_1_len')
parser.add_argument('--Linker_len', default=np.nan)
parser.add_argument('--polyA_2_len', default=np.nan)
parser.add_argument('--read_count')
parser.add_argument('--seq_size')
parser.add_argument('--R1_KnownSeq')
parser.add_argument('--R2_KnownSeq')
parser.add_argument('--read_length', default=301)
parser.add_argument('--window_size_1', default=5)
parser.add_argument('--window_size_2', default=30)
parser.add_argument('--neg_cutoff', default=3)
parser.add_argument('--pos_cutoff', default=3)
parser.add_argument('--round_n', default=5)
parser.add_argument('--sliding_chunkSize', default=10000)
parser.add_argument('--measure_chunkSize', default=500)
parser.add_argument('--max_workers', default=50)

args = parser.parse_args()

sampleID = args.sampleID
input_path = args.input_path
output_path = args.output_path
polyA_1_len = int(args.polyA_1_len)
Linker_len = int(args.Linker_len) if args.Linker_len is not np.nan else np.nan
polyA_2_len = int(args.polyA_2_len) if args.polyA_2_len is not np.nan else np.nan
read_count = int(args.read_count)
seq_size = int(args.seq_size)
R1_KnownSeq = args.R1_KnownSeq
R2_KnownSeq = args.R2_KnownSeq
read_length = int(args.read_length)
window_size_1 = int(args.window_size_1)
window_size_2 = int(args.window_size_2)
neg_cutoff = int(args.neg_cutoff)
pos_cutoff = int(args.pos_cutoff)
round_n = int(args.round_n)
sliding_chunkSize = int(args.sliding_chunkSize)
measure_chunkSize = int(args.measure_chunkSize)
max_workers = int(args.max_workers)

##################################
##### Define Custom Function #####
##################################
###################################################################################
def time_conversion(duration_time):
    hour, remainder = divmod(duration_time, 3600)
    minute, second = divmod(remainder, 60)
    
    return f'{int(hour)}h {int(minute)}m {int(second)}s'
###################################################################################
def define_polyA_structure(sampleID, read_N, read_length, read_count, seq_size, KnownSeq):

    if polyA_2_len is np.nan:
        withLinker, polyA_len_lst = False, [polyA_1_len]
        if read_N == 1: 
            read_region_lst = ["End", "Homopolymer", "UTR", "NA"]
        elif read_N == 2:
            read_region_lst = ["UTR", "Homopolymer", "End", "NA"]
            
    else:
        withLinker, polyA_len_lst = True, [polyA_1_len, Linker_len, polyA_2_len]
        if read_N == 1: 
            polyA_len_lst = list(reversed(polyA_len_lst))
            read_region_lst  = ["End", "Homopolymer_2", "Linker", "Homopolymer_1", "UTR", "NA"]
        
        elif read_N == 2:
            read_region_lst = ["UTR", "Homopolymer_1", "Linker", "Homopolymer_2", "End", "NA"]
    
    #### Make 'vline_lst' 
    read_len_stackedlst, read_len_stackedPos = [len(KnownSeq)], len(KnownSeq)

    for i in range(len(polyA_len_lst)):
        read_len_stackedPos += read_len_stackedlst[i]
        read_len_stackedlst.append(read_len_stackedPos)

    read_len_stackedlst.append(seq_size)
    read_len_stackedlst.append(read_length)


    ### Make 'region_lenDic' 
    read_region_lenDic = {}

    for i in range(len(read_region_lst)):
        read_region_lenDic[read_region_lst[i]] = read_len_stackedlst[i]

    return {'withLinker': withLinker, 
            'polyA_len_lst': polyA_len_lst,
            'polyA_region_lst': read_region_lst[1:-2],
            'read_region_lenDic': read_region_lenDic}
###################################################################################
def Data_bySlidingWindows_chunkProcess(chunk_df, window_size, round_n):
    chunk_columns = chunk_df.columns.tolist()
    
    Q_mean_df = chunk_df[chunk_columns[::-1]].rolling(window=window_size, min_periods=1, axis=1).mean().round(round_n)
    Q_std_df = chunk_df[chunk_columns[::-1]].rolling(window=window_size, min_periods=1, axis=1).std().round(round_n)
    Q_mean_df = Q_mean_df[chunk_columns].astype(np.float32)
    Q_std_df = Q_std_df[chunk_columns].astype(np.float32)
    
    Q_meanDiff_df = Q_mean_df.diff(axis=1).fillna(0).round(round_n).astype(np.float32)
    Q_stdDiff_df = Q_std_df.diff(axis=1).fillna(0).round(round_n).astype(np.float32)
    
    Q_space_meanDiff_df = Q_mean_df.diff(axis=1, periods=window_size).fillna(0).round(round_n).astype(np.float32)
    Q_space_stdDiff_df = Q_std_df.diff(axis=1, periods=window_size).fillna(0).round(round_n).astype(np.float32)
    
    return (Q_meanDiff_df, Q_stdDiff_df, Q_space_meanDiff_df, Q_space_stdDiff_df)

def Data_by_SlidingWindows_chunkedParallel(BaseQ_df, window_size, chunk_size, max_workers, round_n):
    print(f"Start:\t{time.ctime()}"); start_time = time.time()
    
    Q_meanDiff_df_list, Q_stdDiff_df_list = [], []
    Q_space_meanDiff_df_list, Q_space_stdDiff_df_list = [], []

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        num_chunks = (BaseQ_df.shape[0] + chunk_size - 1) // chunk_size
        for start in range(0, BaseQ_df.shape[0], chunk_size):
            end = min(start + chunk_size, BaseQ_df.shape[0])
            chunk_df = BaseQ_df.iloc[start:end]
            futures.append(executor.submit(Data_bySlidingWindows_chunkProcess, chunk_df, window_size, round_n))

        with alive_bar(num_chunks, force_tty=True) as bar:
            for future in as_completed(futures):
                (Q_meanDiff_df, Q_stdDiff_df, Q_space_meanDiff_df, Q_space_stdDiff_df) = future.result()
                Q_meanDiff_df_list.append(Q_meanDiff_df)
                Q_stdDiff_df_list.append(Q_stdDiff_df)
                Q_space_meanDiff_df_list.append(Q_space_meanDiff_df)
                Q_space_stdDiff_df_list.append(Q_space_stdDiff_df)
                bar()

    Q_meanDiff_df = pd.concat(Q_meanDiff_df_list)
    Q_stdDiff_df = pd.concat(Q_stdDiff_df_list)
    Q_space_meanDiff_df = pd.concat(Q_space_meanDiff_df_list)
    Q_space_stdDiff_df = pd.concat(Q_space_stdDiff_df_list)

    print(f"End:\t{time.ctime()}"); end_time = time.time()
    
    duration_time = end_time - start_time; conversion_time = time_conversion(duration_time)
    print(f"Time:\t{duration_time} ({conversion_time})") 
    
    Diff_df_dic = {'Mean': Q_meanDiff_df, 'STdev': Q_stdDiff_df}
    SpaceDiff_df_dic = {'Mean': Q_space_meanDiff_df, 'STdev': Q_space_stdDiff_df}

    return {'Diff': Diff_df_dic, 'SpaceDiff': SpaceDiff_df_dic}
###################################################################################
def Finding_continuousLoc(positions, values, condition, cutoff):
    continuousLoc_lst, collection_lst = [], []
    
    for pos, value in zip(positions, values):
        if condition == 'positive' and value > 0:
            continuousLoc_lst.append(pos)
        elif condition == 'negative' and value < 0:
            continuousLoc_lst.append(pos)
        else:
            if len(continuousLoc_lst) >= cutoff:
                collection_lst.append(continuousLoc_lst)
            continuousLoc_lst = []

    if len(continuousLoc_lst) >= cutoff:
        collection_lst.append(continuousLoc_lst)
    
    return collection_lst

def Define_beforeH_End(window_size, Diff_merge_df, start, collection_cd):
    
    statDF = Diff_merge_df
    
    candidates_1 = [l for lst in collection_cd for l in lst if l >= start]
    
    cd1_statDF = statDF.loc[candidates_1]
    filter1_StatDF = cd1_statDF[(cd1_statDF['STdev'] > 0)]
    
    candidates_2 = [(cd + (window_size - 2)) for cd in filter1_StatDF.index]
    
    if len(candidates_2) > 0:
        return candidates_2[0]
    else:
        return None
    
def Define_Linker_End(window_size, Diff_merge_df, start, collection_cd):
    
    statDF = Diff_merge_df
    
    candidates_1 = [l for lst in collection_cd for l in lst if l >= start]
    
    cd1_statDF = statDF.loc[candidates_1]    
    filter1_StatDF = cd1_statDF[(cd1_statDF['STdev'] < 0)]
    
    candidates_2 = [(cd + (window_size - 1)) for cd in filter1_StatDF.index]

    if len(candidates_2) > 0:
        return candidates_2[0]
    else:
        return None
    
def Define_finalEnd(window_size, Diff_merge_df, SpaceDiff_merge_df, start, end):
    
    statDF = Diff_merge_df
    space_statDF = SpaceDiff_merge_df.loc[start:]
    
    filter1_SpaceStatDF = space_statDF[(space_statDF['Mean'] < 0) & (space_statDF['STdev'] > 0)]
    
    filter2_SpaceStatDF = filter1_SpaceStatDF[filter1_SpaceStatDF['STdev'] >= filter1_SpaceStatDF['STdev'].quantile(0.9)]
    filter2_cd = filter2_SpaceStatDF.index
    filter2_StatDF = statDF.loc[filter2_cd]
    
    filter3_StatDF = filter2_StatDF[filter2_StatDF['Mean'] < 0]
    final_cd = filter3_StatDF.index
    
    if (len(final_cd) == 0) or (final_cd[0] >= end) :
        filter2_SpaceStatDF = filter1_SpaceStatDF[filter1_SpaceStatDF['STdev'] >= filter1_SpaceStatDF['STdev'].quantile(0.75)]
        filter2_cd = filter2_SpaceStatDF.index
        filter2_StatDF = statDF.loc[filter2_cd]

        filter3_StatDF = filter2_StatDF[filter2_StatDF['Mean'] < 0]

        final_cd = filter3_StatDF.index
    
    if (len(final_cd) == 0):
        return None
    
    ampliconEnd = final_cd[0]
    return ampliconEnd

def Measure_PolyA(withLinker, read_idx, read_length, window_size, region_lenDic,
                  w1_Diff_stat_df, w2_SpaceDiff_stat_df,
                  MeanDiff_negLoc_collect, MeanDiff_posLoc_collect):
    
    region_lst = list(region_lenDic.keys())
    
    if withLinker is True:
        beforeH_End = Define_beforeH_End(window_size, w1_Diff_stat_df, region_lenDic[region_lst[0]], MeanDiff_negLoc_collect)
        if beforeH_End is not None:
            beforeH_End_len = beforeH_End - region_lenDic[region_lst[0]]

            Linker_End = Define_Linker_End(window_size, w1_Diff_stat_df, beforeH_End, MeanDiff_posLoc_collect)
            if Linker_End is not None:
                Linker_len = Linker_End - beforeH_End

                afterH_End = Define_finalEnd(window_size, w1_Diff_stat_df, w2_SpaceDiff_stat_df, Linker_End, read_length)
                if afterH_End is not None:
                    afterH_End_len = afterH_End - Linker_End

                else:
                    afterH_End_len = np.nan
            else:
                Linker_len = np.nan
                afterH_End_len = np.nan
        else:
            beforeH_End_len = np.nan
            Linker_len = np.nan
            afterH_End_len = np.nan
            
        merge_len = [beforeH_End_len, Linker_len, afterH_End_len]
    
    else:        
        ampliconEnd = Define_finalEnd(window_size, w1_Diff_stat_df, w2_SpaceDiff_stat_df, region_lenDic[region_lst[0]], read_length)
        if ampliconEnd is not None:
            H_len = ampliconEnd - region_lenDic[region_lst[0]]        
        else:
            H_len = np.nan

        merge_len = [H_len]
    
    return merge_len

def create_chunks_data(data, chunk_size):
    chunks = []
    for i in range(0, len(data), chunk_size):
        df = data.iloc[i:i + chunk_size]
        chunks.append(df)
    return chunks

def create_chunked_idx(idx, chunk_size):
    chunks_idx = []
    for i in range(0, len(idx), chunk_size):
        if (i + chunk_size) > len(idx):
            idx_list = [i for i in range(i, len(idx))]
        else:
            idx_list = [i for i in range(i, i + chunk_size)]
        chunks_idx.append(idx_list)
    return chunks_idx

def Measure_PolyA_chunkProcess(chunk_idx, read_length,
                               withLinker, window_size, region_lenDic,
                               Diff_mean_df_w1, Diff_stdev_df_w1, 
                               SpaceDiff_mean_df_w2, SpaceDiff_stdev_df_w2,
                               MeanDiff_neg_cutoff, MeanDiff_pos_cutoff):
    
    local_lenDic = {}

    for idx in chunk_idx:
        Diff_merge_df_w1 = pd.DataFrame({
            'Mean': Diff_mean_df_w1.loc[idx],
            'STdev': Diff_stdev_df_w1.loc[idx]
        })

        SpaceDiff_merge_df_w2 = pd.DataFrame({
            'Mean': SpaceDiff_mean_df_w2.loc[idx],
            'STdev': SpaceDiff_stdev_df_w2.loc[idx]
        })

        positions = Diff_merge_df_w1.index
        values = Diff_merge_df_w1['Mean'].values

        collect_MeanDiff_negLoc_lst = Finding_continuousLoc(positions, values, 'negative', MeanDiff_neg_cutoff)
        collect_MeanDiff_posLoc_lst = Finding_continuousLoc(positions, values, 'positive', MeanDiff_pos_cutoff)

        polyA_len_list = Measure_PolyA(withLinker, idx, read_length, window_size, region_lenDic,
                                 Diff_merge_df_w1, SpaceDiff_merge_df_w2,
                                 collect_MeanDiff_negLoc_lst, collect_MeanDiff_posLoc_lst)

        local_lenDic[idx] = polyA_len_list

    return local_lenDic

def Measure_PolyA_chunkedParallel(read_idx, BaseQ_w1_dfs_dic, BaseQ_w2_dfs_dic, 
                                  withLinker, read_length, region_lenDic, polyA_region_lst,
                                  window_size, MeanDiff_neg_cutoff, MeanDiff_pos_cutoff,
                                  chunk_size, max_workers, round_n):
    print(f"Start:\t{time.ctime()}"); start_time = time.time()
    
    merge_lenDic = {}
    
    chunked_idx = create_chunked_idx(read_idx, chunk_size)

    Diff_mean_df_w1 = BaseQ_w1_dfs_dic['Diff']['Mean'].sort_index()
    Diff_stdev_df_w1 = BaseQ_w1_dfs_dic['Diff']['STdev'].sort_index()

    SpaceDiff_mean_df_w2 = BaseQ_w2_dfs_dic['SpaceDiff']['Mean'].sort_index()
    SpaceDiff_stdev_df_w2 = BaseQ_w2_dfs_dic['SpaceDiff']['STdev'].sort_index()

    Diff_mean_df_w1_chunks = create_chunks_data(Diff_mean_df_w1, chunk_size)
    Diff_stdev_df_w1_chunks = create_chunks_data(Diff_stdev_df_w1, chunk_size)
    SpaceDiff_mean_df_w2_chunks = create_chunks_data(SpaceDiff_mean_df_w2, chunk_size)
    SpaceDiff_stdev_df_w2_chunks = create_chunks_data(SpaceDiff_stdev_df_w2, chunk_size)
    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        print("Start tasks submission")
        for chunkGroup_n, idx in enumerate(chunked_idx):

            future = executor.submit(Measure_PolyA_chunkProcess, idx, read_length,
                                       withLinker, window_size, region_lenDic,
                                       Diff_mean_df_w1_chunks[chunkGroup_n], Diff_stdev_df_w1_chunks[chunkGroup_n], 
                                       SpaceDiff_mean_df_w2_chunks[chunkGroup_n], SpaceDiff_stdev_df_w2_chunks[chunkGroup_n],
                                       MeanDiff_neg_cutoff, MeanDiff_pos_cutoff)
                
            futures.append(future)
        print("All tasks submitted")

        with alive_bar(len(futures), force_tty=True) as bar:
            for future in as_completed(futures):
                result = future.result()
                merge_lenDic.update(result)
                bar() 

    merge_lenDF = pd.DataFrame.from_dict(merge_lenDic, orient='index').sort_index()
    merge_lenDF.columns = polyA_region_lst
    
    print(f"End:\t{time.ctime()}"); end_time = time.time()
    
    duration_time = end_time - start_time; conversion_time = time_conversion(duration_time)
    print(f"Time:\t{duration_time} ({conversion_time})") 

    return merge_lenDF
###################################################################################

print(f'==================================================================================================================================================================')
print(f"Command :")
print(f"python 2-1_Measurement_polyA_structure.py --prefix {sampleID} --input {input_path} --output {output_path} --polyA_1 {polyA_1_len} --linker {Linker_len} --polyA_2 {polyA_2_len} --readCount {read_count} --seqSize {seq_size} --R1_knownSeq {R1_KnownSeq} --R2_knownSeq {R2_KnownSeq} --readLength {read_length} --windowSize_1 {window_size_1} --windowSize_2 {window_size_2} --neg_cutoff {neg_cutoff} --pos_cutoff {pos_cutoff} --round_n {round_n} --sliding_chunkSize {sliding_chunkSize} --measure_chunkSize {measure_chunkSize} --max_workers {max_workers}")
print(f'==================================================================================================================================================================')
print()

try: 
    if not os.path.exists(f'{output_path}'): 
        os.makedirs(f'{output_path}')
except OSError: 
    print("Error: Failed to create the directory.")
    
print(f'-----------------------------------------------------------------------------------------------------')
print(f"Start:\t{time.ctime()}"); start_time = time.time()
print(f'-----------------------------------------------------------------------------------------------------')


for read_N in [1, 2]:
    if read_N == 1: KnownSeq = R1_KnownSeq
    elif read_N == 2: KnownSeq = R2_KnownSeq
    
    polyA_structure_inf = define_polyA_structure(sampleID, read_N, read_length, read_count, seq_size, KnownSeq)

    print()
    print("######################################################")
    print(f"READ {read_N}")  
    print("######################################################")

    ### 1. Import data
    print("### 1. Importing CSV (Conversion Base Quality Score)")
    BaseQ_df = pd.read_csv(f"{input_path}/{sampleID}_BaseQ_convert.R{read_N}.txt.gz", header=None, sep="\t", dtype='int8', compression='gzip')
    BaseQ_df.columns = [i for i in range(1,BaseQ_df.shape[1]+1)]
    read_idx = list(BaseQ_df.index)
    print()
    
    ### 2. Calculation Statistics Values by Sliding Window Algorithm
    print("### 2. Calculation Statistics Values by Sliding Window Algorithm")
    BaseQ_w1_sliding_data = Data_by_SlidingWindows_chunkedParallel(BaseQ_df, window_size_1, sliding_chunkSize, max_workers, round_n)
    BaseQ_w2_sliding_data = Data_by_SlidingWindows_chunkedParallel(BaseQ_df, window_size_2, sliding_chunkSize, max_workers, round_n)
    print()
    
    ### 3. Measurement Poly(A) Length
    print("### 3.  Measurement Poly(A) Length")
    merge_lenDF = Measure_PolyA_chunkedParallel(read_idx, BaseQ_w1_sliding_data, BaseQ_w2_sliding_data, polyA_structure_inf['withLinker'], read_length, polyA_structure_inf['read_region_lenDic'], polyA_structure_inf['polyA_region_lst'], window_size_1, neg_cutoff, pos_cutoff, measure_chunkSize, max_workers, round_n)
    
    merge_lenDF.to_csv(f'{output_path}/{sampleID}_R{read_N}_polyA_lenDF.csv', sep="\t", encoding='utf-8')
    print("----------------------- Save measured poly(A) length table ------------------------")
    print()
    
print(f'-----------------------------------------------------------------------------------------------------')
print(f"End:\t{time.ctime()}"); end_time = time.time()
print(f'-----------------------------------------------------------------------------------------------------')
print()

duration_time = end_time - start_time
conversion_time = time_conversion(duration_time)
print(f"Time:\t{duration_time} ({conversion_time})")