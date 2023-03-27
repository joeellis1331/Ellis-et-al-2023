#!/usr/bin/python3

import argparse
import traceback
import os
import numpy
import pandas
from statistics import mean
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import math
from scipy import stats
from scipy.stats import ttest_ind
from scipy.optimize import curve_fit
from scipy.stats import f_oneway
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from Bio import SeqIO
import re
from collections import Counter
from maxentpy import maxent
from maxentpy.maxent import load_matrix5, load_matrix3


#set up ArgumentParser
parser = argparse.ArgumentParser(description = "curve fitting")
#section arguments
parser.add_argument("--part1", action = 'store_true', help = "required to run part 1")
parser.add_argument("--part2", action = 'store_true', help = "required to run part 2")
parser.add_argument("--part3", action = 'store_true', help = "required to run part 3")
parser.add_argument("--part4", action = 'store_true', help = "required to run part 4")
parser.add_argument("-a", action = 'store_true')
parser.add_argument("-b", action = 'store_true')
parser.add_argument("-c", action = 'store_true')
parser.add_argument("-d", action = 'store_true')
parser.add_argument("-e", action = 'store_true')
parser.add_argument("-f", action = 'store_true')
#naming shortcut
args = parser.parse_args()


'''
setting the default text to Arial
'''
matplotlib.rcParams['font.family'] = 'sans-serif'


'''this section reads in dataframes, the idea is to first filter the d1000 dataset to get the largest set of events because MBNL
    levels are highest, and then simply go through and find the matching events in the other data. Save final dataset for use later'''
if args.part1 == True:
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

    '''does all the appropriate filtering for the SE.JCEC.txt files'''
    def filter(file, suf_con, suf_dose): #file, suffix control, suffix dose
        df = pandas.read_table(file) #create input df

        col_list = ['IJC_SAMPLE_1', 'SJC_SAMPLE_1', 'IJC_SAMPLE_2', 'SJC_SAMPLE_2', 'IncLevel1', 'IncLevel2']

        for col in col_list:
            df[col] = df.apply(lambda x: convert(x[col]), axis = 1) #converts columns with triplicates to float list
            df['%s_avg' % col] = df.apply(lambda x: round(mean(x[col]), 2), axis = 1) #averages IncLevel

        #only filters df1000 dataset, this is the starting events since it has the largest change of diff splicing
        if '1000' in suf_dose:
            df = df[(df.FDR < 0.05) & (df.IncLevelDifference.abs() >= 0.2)] #input filtering (FDR, PSI, and read counts)
            df = df[((df.IJC_SAMPLE_1_avg + df.SJC_SAMPLE_1_avg) >= 5.0) & ((df.IJC_SAMPLE_2_avg + df.SJC_SAMPLE_2_avg) >= 5.0)] #filters sum of reads
        else:
            pass

        #subset data and remove extraneous columns
        df = df[['GeneID', 'geneSymbol', 'chr', 'strand', 'exonStart_0base',
                'exonEnd', 'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE', 'IJC_SAMPLE_1_avg', 'SJC_SAMPLE_1_avg', 'IJC_SAMPLE_2_avg', 'SJC_SAMPLE_2_avg',
                'IncLevel1', 'IncLevel1_avg', 'IncLevel2', 'IncLevel2_avg']]
        df = df.rename(columns = {'IncLevel1_avg':('IncLevel_avg' + suf_con), 'IncLevel1':('IncLevel' + suf_con), 'IncLevel2_avg':('IncLevel_avg' + suf_dose), 'IncLevel2':('IncLevel' + suf_dose),
                                 'IJC_SAMPLE_1_avg':('IJC_SAMPLE_1_avg' + suf_con), 'SJC_SAMPLE_1_avg':('SJC_SAMPLE_1_avg' + suf_con), 'IJC_SAMPLE_2_avg':('IJC_SAMPLE_2_avg' + suf_dose), 'SJC_SAMPLE_2_avg':('SJC_SAMPLE_2_avg' + suf_dose)})

        return df

    path = '/mnt/c/Users/joeel/Desktop/Transient_Coreg_Files/NEW_bioinfo/rMATS_Results/new_ddHEK_triplicate/'
    df100 = filter(path + 'd0_vs_d100_HEK/SE.MATS.JCEC.txt', '_d0', '_d100')
    df170 = filter(path + 'd0_vs_d170_HEK/SE.MATS.JCEC.txt', '_d0', '_d170')
    df250 = filter(path + 'd0_vs_d250_HEK/SE.MATS.JCEC.txt', '_d0', '_d250')
    df280 = filter(path + 'd0_vs_d280_HEK/SE.MATS.JCEC.txt', '_d0', '_d280')
    df1000 = filter(path + 'd0_vs_d1000_HEK/SE.MATS.JCEC.txt', '_d0', '_d1000')

    #merge all individual pairwise comparisons into one big dataframe
    merge_list = ['GeneID', 'geneSymbol', 'chr', 'strand', 'exonStart_0base',
           'exonEnd', 'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE', 'IncLevel_avg_d0']
    df_all = df1000.merge(df100, how = 'outer', on = merge_list).dropna()
    df_all = df_all.merge(df170, how = 'outer', on = merge_list).dropna()
    df_all = df_all.merge(df250, how = 'outer', on = merge_list).dropna()
    df_all = df_all.merge(df280, how = 'outer', on = merge_list).dropna()
    df_all = df_all.drop_duplicates(subset = ['geneSymbol', 'exonStart_0base', 'exonEnd'])

    #reorganize columns
    df_all = df_all[['GeneID', 'geneSymbol', 'chr', 'strand', 'exonStart_0base', 'exonEnd', 'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE',
        'IJC_SAMPLE_1_avg_d0', 'SJC_SAMPLE_1_avg_d0', 'IJC_SAMPLE_2_avg_d100', 'SJC_SAMPLE_2_avg_d100', 'IJC_SAMPLE_2_avg_d170', 'SJC_SAMPLE_2_avg_d170',
        'IJC_SAMPLE_2_avg_d250', 'SJC_SAMPLE_2_avg_d250', 'IJC_SAMPLE_2_avg_d280', 'SJC_SAMPLE_2_avg_d280',
        'IJC_SAMPLE_2_avg_d1000', 'SJC_SAMPLE_2_avg_d1000', 'IncLevel_d0', 'IncLevel_d100', 'IncLevel_d170', 'IncLevel_d250', 'IncLevel_d280', 'IncLevel_d1000',
        'IncLevel_avg_d0', 'IncLevel_avg_d100', 'IncLevel_avg_d170', 'IncLevel_avg_d250', 'IncLevel_avg_d280', 'IncLevel_avg_d1000']]

    #save as csv and pickle for viewing or later use
    df_all.to_csv('MBNLeventlist_original_HEK.csv', index = False)
    df_all.to_pickle('MBNLeventlist_original_HEK.pkl')

'''curve fit section'''
if args.part2 == True:
    '''curve fits and associated parameters'''
    if args.a == True:
        df_all = pandas.read_pickle('/mnt/c/Users/joeel/Desktop/Transient_Coreg_Files/NEW_bioinfo/HEK/MBNLeventlist_original_HEK.pkl')

        #change IncLevel avg colmumn values to 0-100 scale instead of 0-1
        def percentify(list):
            new_list = []
            for i in list:
                i = round(i * 100, 1)
                new_list.append(i)
            return new_list

        df_all['IncLevel_d0'] = df_all.apply(lambda x: percentify(x['IncLevel_d0']), axis = 1)
        df_all['IncLevel_d100'] = df_all.apply(lambda x: percentify(x['IncLevel_d100']), axis = 1)
        df_all['IncLevel_d170'] = df_all.apply(lambda x: percentify(x['IncLevel_d170']), axis = 1)
        df_all['IncLevel_d250'] = df_all.apply(lambda x: percentify(x['IncLevel_d250']), axis = 1)
        df_all['IncLevel_d280'] = df_all.apply(lambda x: percentify(x['IncLevel_d280']), axis = 1)
        df_all['IncLevel_d1000'] = df_all.apply(lambda x: percentify(x['IncLevel_d1000']), axis = 1)

        df_all['IncLevel_avg_d0'] = round(df_all['IncLevel_avg_d0'] * 100, 1)
        df_all['IncLevel_avg_d100'] =  round(df_all['IncLevel_avg_d100'] * 100, 1)
        df_all['IncLevel_avg_d170'] = round(df_all['IncLevel_avg_d170'] * 100, 1)
        df_all['IncLevel_avg_d250'] = round(df_all['IncLevel_avg_d250'] * 100, 1)
        df_all['IncLevel_avg_d280'] = round(df_all['IncLevel_avg_d280'] * 100, 1)
        df_all['IncLevel_avg_d1000'] = round(df_all['IncLevel_avg_d1000'] * 100, 1)

        #sigmoid curve fit function
        '''NOTE: The LogEC50 parameter matches the x axis no transformation needed'''
        def sigmoid(x, top, bottom, LogEC50, slope):
            y = bottom + ((top - bottom)/(1 + 10**((LogEC50 - x) * slope)))
            return (y)

        def get_curve_data(row):
            '''NOTE:
                -x data is log(dose + 1) so the 0 point can be used
                -the final cellular concentration is 100 fold less (ie. d100 --> d1)'''
            #set x and y data
            xdata = [math.log10(1), math.log10(2), math.log10(2.7), math.log10(3.5), math.log10(3.8), math.log10(11)]
            ydata = [row.IncLevel_avg_d0, row.IncLevel_avg_d100, row.IncLevel_avg_d170, row.IncLevel_avg_d250, row.IncLevel_avg_d280, row.IncLevel_avg_d1000]

            ############## find standard deviation for each ydata point #####################################
            yerr = [numpy.std(row.IncLevel_d0), numpy.std(row.IncLevel_d100), numpy.std(row.IncLevel_d170),
                    numpy.std(row.IncLevel_d250), numpy.std(row.IncLevel_d280), numpy.std(row.IncLevel_d1000)]

            #calculate pooled StDev
            sq_sum_std = 0
            #takes the square of each std and sums them
            for std in yerr:
                sq_sum_std += std**2
            #the actual pooled std calculation
            pool = round((math.sqrt(sq_sum_std/2)), 3)

            ############# dose response model  curve fit ############################
            try:
                #p0 is initial guesses
                if row.IncLevel_avg_d1000 > row.IncLevel_avg_d0:
                    p0 = [max(ydata), min(ydata), 0.5, 3] #mean of log dox values is 0.48, slope mean from prev boxplots is ~3
                else:
                    p0 = [min(ydata), max(ydata), 0.5, 3] #mean of log dox values is 0.48, slope mean from prev boxplots is ~3
                #fit curve
                popt, pcov = curve_fit(sigmoid, xdata, ydata, p0)

            except Exception:
                #ignores inability to fit a curve error but prints other errors
                error = traceback.format_exc()
                if 'Optimal parameters not found' in error:
                    pass
                else:
                    print(error)

            else:
                #residual sum of squares calculations
                residuals = ydata - sigmoid(xdata, *popt)
                ss_res = numpy.sum(residuals**2)
                #total sum of squares
                ss_tot = numpy.sum((ydata-numpy.mean(ydata))**2)
                #r sq value
                r_squared = 1 - (ss_res/ss_tot)
                #store these values in the dataframe
                df_all.at[row.Index, 'RSS'] = ss_res
                df_all.at[row.Index, 'total_sumsq'] = ss_res
                df_all.at[row.Index, 'R_sq'] = r_squared
                df_all.at[row.Index, 'LogRSS'] = math.log10(ss_res)

                #standard error calculation, list (standard deviation / sqrt(# samples))
                stderr = numpy.sqrt(numpy.diag(pcov))/numpy.sqrt(len(ydata))
                #appends paramters to DataFrame
                df_all.at[row.Index, 'stderr_top'] = round(stderr[0], 1)
                df_all.at[row.Index, 'stderr_bottom'] = round(stderr[1], 1)
                df_all.at[row.Index, 'stderr_LogEC50'] = round(stderr[2], 2)
                df_all.at[row.Index, 'stderr_slope'] = round(stderr[3], 1)

                #appends paramters to DataFrame
                df_all.at[row.Index, 'top'] = max([popt[0], popt[1]])
                df_all.at[row.Index, 'bottom'] = min([popt[0], popt[1]])
                df_all.at[row.Index, 'LogEC50'] = popt[2]
                df_all.at[row.Index, 'slope'] = popt[3]
                df_all.at[row.Index, 'pooled_StDev'] = pool

        #################################### running the actual function #####################################
        for row in df_all.itertuples():
            get_curve_data(row)

        #this removes any sections that didnt get a curve fit because the slope/EC50 weren't added but instead NAs are present
        df_all = df_all.dropna()
        #final fitstats dataframe
        df_all.to_pickle('MBNLeventlist_curvefit_HEK.pkl')

    '''plot the actual curves based on filtering criteria'''
    if args.b == True:
        #sigmoid curve fit function
        '''note! The LogEC50 parameter matches the x axis no transformation needed. for whatever reason calculating the log
        as part of the equation messes up the fitting stats'''
        def sigmoid(x, top, bottom, LogEC50, slope):
            y = bottom + ((top - bottom)/(1 + 10**((LogEC50 - x) * slope)))
            return (y)
    
        #takes parameters from text file and filters based on something, then plots curves
        def plot_curves(ind):
            xdata = [math.log10(1), math.log10(2), math.log10(2.7), math.log10(3.5), math.log10(3.8), math.log10(11)]
            ydata = [df_all.iloc[ind].IncLevel_avg_d0, df_all.iloc[ind].IncLevel_avg_d100, df_all.iloc[ind].IncLevel_avg_d170, df_all.iloc[ind].IncLevel_avg_d250, df_all.iloc[ind].IncLevel_avg_d280, df_all.iloc[ind].IncLevel_avg_d1000]

            ############## find std dev #####################################
            yerr = [numpy.std(df_all.iloc[ind].IncLevel_d0), numpy.std(df_all.iloc[ind].IncLevel_d100), numpy.std(df_all.iloc[ind].IncLevel_d170),
                    numpy.std(df_all.iloc[ind].IncLevel_d250), numpy.std(df_all.iloc[ind].IncLevel_d280), numpy.std(df_all.iloc[ind].IncLevel_d1000)]

            #define plus minus symbol for plotting
            plus_minus = (u'\u00B1')

            #RNAseq 
            #use 3 filters before plotting actual curves
            #residual sum of squares cutoff
            if df_all.iloc[ind].LogRSS <= 1.5:
                #cutoff for top and bottom to be realistic
                if (df_all.iloc[ind].top <= 100.0) and (df_all.iloc[ind].bottom >= 0.0):
                    #cutoff for slope, lets say 16, thats 2 times stacy's published RT for MBNL2
                    if df_all.iloc[ind].slope < 16.0:
                        #plots actual datapoints with error bars
                        plt.plot(xdata, ydata, 'o', c = 'black')
                        plt.errorbar(xdata, ydata, yerr = yerr, ls = 'none', c = 'black', alpha = 0.5, capsize = 4) #plots errors bars

                        #plot the curve
                        x = numpy.linspace(-0.5, 1.5, num = 100) #returns evenly spaced numbers over a certain interval (start, stop, num = int)
                        if df_all.iloc[ind].IncLevel_avg_d1000 > df_all.iloc[ind].IncLevel_avg_d0:
                            popt = [df_all.iloc[ind].top, df_all.iloc[ind].bottom, df_all.iloc[ind].LogEC50, df_all.iloc[ind].slope]
                        else:
                            popt = [df_all.iloc[ind].bottom, df_all.iloc[ind].top, df_all.iloc[ind].LogEC50, df_all.iloc[ind].slope]
                        y = sigmoid(x, *popt) #fits sigmoid curve to y data
                        plt.plot(x, y,  c = 'orange', label = 'RNAseq fit') #plots sigmoid fit

                        #appends parameters for table plotting to array
                        top_table = f'{round(max([popt[0], popt[1]]), 1)}{plus_minus}{df_all.iloc[ind].stderr_top}'
                        bottom_table = f'{round(min([popt[0], popt[1]]), 1)}{plus_minus}{df_all.iloc[ind].stderr_bottom}'
                        ec50_table = f'{round(popt[3], 2)}{plus_minus}{df_all.iloc[ind].stderr_LogEC50}'
                        slope_table = f'{round(popt[2], 2)}{plus_minus}{df_all.iloc[ind].stderr_slope}'
                        r_sq_table = f'{round(df_all.iloc[ind].R_sq, 3)}'
                        cellText = [[top_table, bottom_table, ec50_table, slope_table, r_sq_table]] #need extra brackets to create 2D array
                        
                        #plots table
                        colLabels = ['$\\bf{Top}$', '$\\bf{Bottom}$', '$\\bf{Slope}$', '$\\bf{LogEC50}$', '$\\bf{R\u00b2}$'] #uses TeX code
                        colColours = ['lightgrey' for x in range(len(colLabels))]
                        tabl = plt.table(cellText = cellText, cellLoc = 'center', colLabels = colLabels, colColours = colColours, colLoc = 'center', loc = 'top')
                        #allows me to set fontsize
                        tabl.auto_set_font_size(False)
                        tabl.set_fontsize(11.5)
                        #changing cell heights
                        tabl_props = tabl.properties()
                        tabl_cells = tabl_props['children']
                        for cell in tabl_cells: cell.set_height(0.075)

                        #final plot adjustments
                        plt.subplots_adjust(top = 0.80)
                        plt.legend(loc = 'best')
                        plt.ylabel('PSI', fontsize = 'xx-large')
                        plt.xlabel('log[(dox + 1)]', fontsize = 'xx-large')
                        plt.xticks(ticks = [-0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5], labels = [-0.5, '', 0, '', 0.5, '', 1.0, '', 1.5], fontsize = 17)
                        plt.yticks(fontsize = 17)
                        plt.title(f'{df_all.iloc[ind].geneSymbol}', pad = 50, fontsize = 'xx-large', style = 'italic')
                        plt.suptitle(f'exon coord:{df_all.iloc[ind].exonStart_0base}-{df_all.iloc[ind].exonEnd}', x = 0.55, y = 0.84, fontsize = 'xx-large')
                        plt.tight_layout()
                        name = f'{df_all.iloc[ind].geneSymbol}_{df_all.iloc[ind].exonStart_0base}-{df_all.iloc[ind].exonEnd}'
                        plt.savefig(f'./plots/curve_fits/{name}.tiff', dpi = 300, format = 'tiff')
                        plt.close()

        
        #reads in dataframe
        df_all = pandas.read_pickle('/mnt/c/Users/joeel/Desktop/Transient_Coreg_Files/NEW_bioinfo/HEK/MBNLeventlist_curvefit_HEK.pkl')
        #actual function call
        ind_list = [x for x in range(len(df_all))] #prob a better way to do this
        for ind in ind_list:
            plot_curves(ind)
