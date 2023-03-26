#!/usr/bin/python3

import argparse
import numpy
import pandas
import statistics
import seaborn
import matplotlib
import matplotlib.pyplot as plt

#set up ArgumentParser
parser = argparse.ArgumentParser(description = "AS general overview")
#section arguments
parser.add_argument("-a", action = 'store_true')
parser.add_argument("-b", action = 'store_true')
parser.add_argument("-c", action = 'store_true')
#naming shortcut
args = parser.parse_args()


'''
takes each rMATS.JCEC file for all the AS types across the 5 pairwise comparisons, filters them based on reads and IncLevelDiff, then
combines into one dataframe and stores as a pickle for later use
'''
if args.a == True:
    '''converts original rMATS output into a usable list of floats'''
    def convert(string):
            _list = list(string.split(',')) #splits string on comma into elements and stores in list
            while True:
                if 'NA' in _list:
                    _list.remove('NA') #removes any NA's
                else:
                    break
            _list = list(map(float, _list)) #converts str into floats
            return _list

    '''Filters each events by reads and PSI value'''
    def filter_AS(file):
        df = pandas.read_table(file)

        col_list = ['IJC_SAMPLE_1', 'SJC_SAMPLE_1', 'IJC_SAMPLE_2', 'SJC_SAMPLE_2', 'IncLevel1', 'IncLevel2']
        for col in col_list:
            df[col] = df.apply(lambda x: convert(x[col]), axis = 1) #converts columns with triplicates to float list
            df['%s_avg' % col] = df.apply(lambda x: round(statistics.mean(x[col]), 2), axis = 1) #averages IncLevel and read counts for samples

        df = df[(df.FDR <= 0.05) & (df.IncLevelDifference.abs() >= 0.1)] #input filtering (FDR, PSI, and read counts)
        df = df[((df.IJC_SAMPLE_1_avg + df.SJC_SAMPLE_1_avg) >= 5.0) & ((df.IJC_SAMPLE_2_avg + df.SJC_SAMPLE_2_avg) >= 5.0)] #filters sum of reads
        
        #add column with specific type of AS event to sort later
        if 'SE' in file:
            df['AS_type'] = 'SE'
        elif 'MXE' in file:
            df['AS_type'] = 'MXE'
        elif 'A5SS' in file:
            df['AS_type'] = 'A5SS' 
        elif 'A3SS' in file:
            df['AS_type'] = 'A3SS'
        else:
            df['AS_type'] = 'RI'    
        
        return df

    '''collects the dfs for each type of AS event and turns into one big dataframe'''
    def organize_data(dir_name):
        path = '/mnt/c/Users/joeel/Desktop/Transient_Coreg_Files/NEW_bioinfo/rMATS_Results/new_ddHEK_triplicate/'
        abs_path = path + dir_name + '/'

        df_se = filter_AS((abs_path + 'SE.MATS.JCEC.txt'))
        df_mxe = filter_AS((abs_path + 'MXE.MATS.JCEC.txt'))
        df_a5ss = filter_AS((abs_path + 'A5SS.MATS.JCEC.txt'))
        df_a3ss = filter_AS((abs_path + 'A3SS.MATS.JCEC.txt'))
        df_ri = filter_AS((abs_path + 'RI.MATS.JCEC.txt'))

        df = pandas.concat([df_se, df_mxe, df_a3ss, df_a5ss, df_ri])
        return df

    #calls functions and saves dataframes
    df_d100 = organize_data('d0_vs_d100_HEK')
    df_d170 = organize_data('d0_vs_d170_HEK')
    df_d250 = organize_data('d0_vs_d250_HEK')
    df_d280 = organize_data('d0_vs_d280_HEK')
    df_d1000 = organize_data('d0_vs_d1000_HEK')

    #further adds another column noting the dose
    df_d100['dose'], df_d170['dose'], df_d250['dose'], df_d280['dose'], df_d1000['dose'] = 'd100', 'd170', 'd250', 'd280', 'd1000'
    df = pandas.concat([df_d100, df_d170, df_d250, df_d280, df_d1000])
    df.to_pickle('ASevents_all.pkl')

'''make some plots to assess overall AS outcomes'''
if args.b == True:
    #setting the default text to Arial
    matplotlib.rcParams['font.family'] = 'sans-serif'
    #less verbose way to save figs
    plot_path = '/mnt/c/Users/joeel/Desktop/Transient_Coreg_Files/NEW_bioinfo/HEK/plots/ASevent_overview'
    #load data
    df = pandas.read_pickle('/mnt/c/Users/joeel/Desktop/Transient_Coreg_Files/NEW_bioinfo/HEK/ASevents_all.pkl')

    #counts plot
    seaborn.countplot(data = df, x = 'dose', hue = 'AS_type')
    plt.xlabel('d0 vs ___', size = 'xx-large')
    plt.ylabel('Total Events', size = 'xx-large')
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.tight_layout()
    plt.savefig(f'{plot_path}/AS_totals.png', dpi = 600, format = 'png')
    plt.close()

    


