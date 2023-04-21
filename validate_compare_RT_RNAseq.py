#!/usr/bin/python3

import argparse
import numpy
import pandas
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import math
import statistics
import scipy
from scipy.optimize import curve_fit
from Bio import SeqIO
from collections import Counter
from maxentpy import maxent
from maxentpy.maxent import load_matrix5, load_matrix3
import itertools
from scipy.stats import norm

#set up ArgumentParser
parser = argparse.ArgumentParser(description = "post curve fit stuff")
#section arguments
parser.add_argument("--part1", action = 'store_true', help = "required to run part 1")
parser.add_argument("--part2", action = 'store_true', help = "required to run part 2")
parser.add_argument("--part3", action = 'store_true', help = "required to run part 3")
parser.add_argument("--part4", action = 'store_true', help = "required to run part 4")
parser.add_argument("--part5", action = 'store_true', help = "required to run part 5")
parser.add_argument("-a", action = 'store_true')
parser.add_argument("-b", action = 'store_true')
parser.add_argument("-c", action = 'store_true')
parser.add_argument("-d", action = 'store_true')
#naming shortcut
args = parser.parse_args()

'''
setting the default text to Arial
'''
matplotlib.rcParams['font.family'] = 'sans-serif'

'''Generates curve fits first and saves to file
uses that data to then plot curves for original events (no dilutions etc) with table of values on top

The structure of the input excel file (df_all) is paired rows, where the first is the RNAseq and the second is the RT-PCR data for a single event.
'''
if args.part1 == True:
    df_all = pandas.read_excel('/mnt/c/Users/joeel/Desktop/Transient_Coreg_Files/NEW_bioinfo/HEK/MBNLeventlist_curveValidation_HEK_original.xlsx', sheet_name = 0, index_col = 0)
    df_all = df_all.dropna()

    #needs wonky conversion from xlsx file to run in this script
    def convert_to_list(s):
        s = s[1:-1]
        s = s.split(',')
        s = [float(x) for x in s]
        return s

    cols = ['IncLevel_d0', 'IncLevel_d100', 'IncLevel_d170',  'IncLevel_d250', 'IncLevel_d280', 'IncLevel_d1000']
    for col in cols:
        df_all[col] = df_all.apply(lambda x: convert_to_list(x[col]), axis = 1)

    #sigmoid curve fit function
    '''note! The LogEC50 parameter matches the x axis no transformation needed. for whatever reason calculating the log
    as part of the equation messes up the fitting stats'''
    def sigmoid(x, top, bottom, LogEC50, slope):
        y = bottom + ((top - bottom)/(1 + 10**((LogEC50 - x) * slope)))
        return (y)

    #generates curve fits and saves parameters into a dataframe
    def get_curve_data(row):
        '''x data is log(dose + 1) so the 0 point can be used
        also the final cellular concentration is 100 fold less! couldnt name variables with a decimal! (eg. d100 --> d1)'''
        #set x and y data
        xdata = list(itertools.repeat(math.log10(1), len(row.IncLevel_d0))) + list(itertools.repeat(math.log10(2), len(row.IncLevel_d100)))\
            + list(itertools.repeat(math.log10(2.7), len(row.IncLevel_d170))) + list(itertools.repeat(math.log10(3.5), len(row.IncLevel_d250)))\
            + list(itertools.repeat(math.log10(3.8), len(row.IncLevel_d280))) + list(itertools.repeat(math.log10(11), len(row.IncLevel_d1000)))
        ydata = row.IncLevel_d0 + row.IncLevel_d100 + row.IncLevel_d170 + row.IncLevel_d250 + row.IncLevel_d280 + row.IncLevel_d1000
       
        ############## find std dev for each set of samples #####################################
        yerr = [numpy.std(row.IncLevel_d0), numpy.std(row.IncLevel_d100), numpy.std(row.IncLevel_d170),
                numpy.std(row.IncLevel_d250), numpy.std(row.IncLevel_d280), numpy.std(row.IncLevel_d1000)]

        #calculate pooled StDev
        sq_sum_std = 0
        #takes the square of each std and sums them
        for std in yerr:
            sq_sum_std += std**2
        #the actual pooled std calculation
        pool = round((math.sqrt(sq_sum_std/2)), 3)

        ############# dose response model fit ############################
        #the curve fit function doesn't work for all data
        #so this allows me to ignore any errors risen from inability to generate a curve
        #it only plots if a curve is successfully fit
        try:
            #p0 is initial guesses
            if row.IncLevel_avg_d1000 > row.IncLevel_avg_d0:
                p0 = [max(ydata), min(ydata), 0.5, 3] #mean of log dox values is 0.48, slope mean from prev boxplots is like 3
            else:
                p0 = [min(ydata), max(ydata), 0.5, 3] #mean of log dox values is 0.48, slope mean from prev boxplots is like 3

            #fit curve
            popt, pcov = curve_fit(sigmoid, xdata, ydata, p0)

        except BaseException as error:
            print('An exception occurred: {}'.format(error))

        else:
            #general goodness of fit calculations
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

             #confidence interval calculations
            def calc_ci(ci, param, stderr):
                z = norm.ppf((1.0 + ci)/2)
                lower_ci = param - (z * stderr)
                upper_ci = param + (z * stderr)
                return lower_ci, upper_ci
            
            df_all.at[row.Index, 'top_lowCI'], df_all.at[row.Index, 'top_upCI'] = calc_ci(0.95, (max([popt[0], popt[1]])), stderr[0])
            df_all.at[row.Index, 'bottom_lowCI'], df_all.at[row.Index, 'bottom_upCI'] = calc_ci(0.95, (min([popt[0], popt[1]])), stderr[1])
            df_all.at[row.Index, 'LogEC50_lowCI'], df_all.at[row.Index, 'LogEC50_upCI'] = calc_ci(0.95, popt[2], stderr[2])
            df_all.at[row.Index, 'slope_lowCI'], df_all.at[row.Index, 'slope_upCI'] = calc_ci(0.95, popt[3], stderr[3])

            #appends paramters to DataFrame
            df_all.at[row.Index, 'top'] = max([popt[0], popt[1]])
            df_all.at[row.Index, 'bottom'] = min([popt[0], popt[1]])
            df_all.at[row.Index, 'LogEC50'] = popt[2]
            df_all.at[row.Index, 'slope'] = popt[3]
            df_all.at[row.Index, 'pooled_StDev'] = pool
    
    #runs the above function
    for row in df_all.itertuples():
        get_curve_data(row)

    ######################################################################################
    #makes plots
    #####################################################################################
    def plot_curves(ind, cellText):
        xdata = [math.log10(1), math.log10(2), math.log10(2.7), math.log10(3.5), math.log10(3.8), math.log10(11)]
        ydata = [df_all.iloc[ind].IncLevel_avg_d0, df_all.iloc[ind].IncLevel_avg_d100, df_all.iloc[ind].IncLevel_avg_d170, df_all.iloc[ind].IncLevel_avg_d250, df_all.iloc[ind].IncLevel_avg_d280, df_all.iloc[ind].IncLevel_avg_d1000]

        ############## find std dev #####################################
        yerr = [numpy.std(df_all.iloc[ind].IncLevel_d0), numpy.std(df_all.iloc[ind].IncLevel_d100), numpy.std(df_all.iloc[ind].IncLevel_d170),
                numpy.std(df_all.iloc[ind].IncLevel_d250), numpy.std(df_all.iloc[ind].IncLevel_d280), numpy.std(df_all.iloc[ind].IncLevel_d1000)]

        #RNAseq data
        if df_all.iloc[ind].method == 'RNAseq':
            #plots actual datapoints with error bars
            plt.plot(xdata, ydata, 'o', c = 'black', fillstyle = 'none')
            plt.errorbar(xdata, ydata, yerr = yerr, ls = 'none', c = 'darkorange', alpha = 0.75, capsize = 5) #plots errors bars

            #plot the curve
            x = numpy.linspace(-0.5, 1.5, num = 100) #returns evenly spaced numbers over a certain interval (start, stop, num = int)
            if df_all.iloc[ind].IncLevel_avg_d1000 > df_all.iloc[ind].IncLevel_avg_d0:
                popt = [df_all.iloc[ind].top, df_all.iloc[ind].bottom, df_all.iloc[ind].LogEC50, df_all.iloc[ind].slope]
            else:
                popt = [df_all.iloc[ind].bottom, df_all.iloc[ind].top, df_all.iloc[ind].LogEC50, df_all.iloc[ind].slope]
            y = sigmoid(x, *popt) #fits sigmoid curve to y data
            plt.plot(x, y,  c = 'darkorange', label = '%s fit' % df_all.iloc[ind].method) #plots sigmoid fit

            #appends parameters for table plotting to array
            top_table = f'{round(max([popt[0], popt[1]]), 1)}'
            bottom_table = f'{round(min([popt[0], popt[1]]), 1)}'
            ec50_table = f'{round(popt[3], 2)}'
            slope_table = f'{round(popt[2], 1)}'
            r_sq_table = f'{round(df_all.iloc[ind].R_sq, 2)}'
            top_ci = f'{round(df_all.iloc[ind].top_lowCI, 1)}-{round(df_all.iloc[ind].top_upCI, 1)}'
            bottom_ci = f'{round(df_all.iloc[ind].bottom_lowCI, 1)}-{round(df_all.iloc[ind].bottom_upCI, 1)}'
            ec50_ci = f'{round(df_all.iloc[ind].LogEC50_lowCI, 2)}-{round(df_all.iloc[ind].LogEC50_upCI, 2)}'
            slope_ci = f'{round(df_all.iloc[ind].slope_lowCI, 1)}-{round(df_all.iloc[ind].slope_upCI, 1)}'
            extend_list = ['RNAseq', 'Best-fit', top_table, bottom_table, ec50_table, slope_table, r_sq_table]
            cellText[0].extend(extend_list)
            extend_list = [None, '95% CI', top_ci, bottom_ci, slope_ci, ec50_ci, None]
            cellText[1].extend(extend_list)

        #RT-PCR data
        else:
            #plots actual datapoints with error bars
            plt.plot(xdata, ydata, 'o', c = 'black', fillstyle = 'none')
            plt.errorbar(xdata, ydata, yerr = yerr, ls = 'none', c = 'purple', alpha = 0.75, capsize = 5) #plots errors bars

            #plot the curve
            x = numpy.linspace(-0.5, 1.5, num = 100) #returns evenly spaced numbers over a certain interval (start, stop, num = int)
            if df_all.iloc[ind].IncLevel_avg_d1000 > df_all.iloc[ind].IncLevel_avg_d0:
                popt = [df_all.iloc[ind].top, df_all.iloc[ind].bottom, df_all.iloc[ind].LogEC50, df_all.iloc[ind].slope]
            else:
                popt = [df_all.iloc[ind].bottom, df_all.iloc[ind].top, df_all.iloc[ind].LogEC50, df_all.iloc[ind].slope]
            y = sigmoid(x, *popt) #fits sigmoid curve to y data
            plt.plot(x, y,  c = 'purple', label = '%s fit' % df_all.iloc[ind].method) #plots sigmoid fit

            #appends parameters for table plotting to array
            top_table = f'{round(max([popt[0], popt[1]]), 1)}'
            bottom_table = f'{round(min([popt[0], popt[1]]), 1)}'
            ec50_table = f'{round(popt[3], 2)}'
            slope_table = f'{round(popt[2], 1)}'
            r_sq_table = f'{round(df_all.iloc[ind].R_sq, 2)}'
            top_ci = f'{round(df_all.iloc[ind].top_lowCI, 1)}-{round(df_all.iloc[ind].top_upCI, 1)}'
            bottom_ci = f'{round(df_all.iloc[ind].bottom_lowCI, 1)}-{round(df_all.iloc[ind].bottom_upCI, 1)}'
            ec50_ci = f'{round(df_all.iloc[ind].LogEC50_lowCI, 2)}-{round(df_all.iloc[ind].LogEC50_upCI, 2)}'
            slope_ci = f'{round(df_all.iloc[ind].slope_lowCI, 1)}-{round(df_all.iloc[ind].slope_upCI, 1)}'
            extend_list = ['RT-PCR', 'Best-fit', top_table, bottom_table, ec50_table, slope_table, r_sq_table]
            cellText[2].extend(extend_list)
            extend_list = [None, '95% CI', top_ci, bottom_ci, slope_ci, ec50_ci, None]
            cellText[3].extend(extend_list)

        name = '%s_%s-%s' % (df_all.iloc[ind].geneSymbol, str(df_all.iloc[ind].exonStart_0base), str(df_all.iloc[ind].exonEnd))
        return name, cellText


    #necessary plotting stuff, goes by tuple index
    ind = [x for x in range(len(df_all))] #creates list of indexes for len of dataframe
    it = iter(ind) #creates iterator
    ind = list(zip(it, it)) #accesses iterator and advances position, creating a tuple of pairs e.g. ((1,2), (3,4))
   
    #iterates through index tuple
    for tup in ind:
        #2D list for plotting table
        #cellText = [['RNAseq'], ['RT-PCR']]
        cellText = [[], [], [], []]
        #plots and grabs geneSymbol for chart title
        name, cellText = plot_curves(tup[0], cellText) #RNAseq
        name, cellText = plot_curves(tup[1], cellText) #RT-PCR

        #plots table
        colLabels = ['$\\bf{Method}$', None, '$\\bf{Top}$', '$\\bf{Bottom}$', '$\\bf{Slope}$', '$\\bf{LogEC50}$', '$\\bf{R\u00b2}$'] #uses TeX code
        colColours = ['lightgrey' for x in range(len(colLabels))]
        tabl = plt.table(cellText = cellText, cellLoc = 'center', colLabels = colLabels, colColours = colColours, colLoc = 'center', loc = 'top')
        #allows me to set fontsize
        tabl.auto_set_font_size(False)
        tabl.set_fontsize(10.5)
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
        plt.title(df_all.iloc[tup[0]].geneSymbol, pad = 75, fontsize = 'xx-large', style = 'italic')
        plt.tight_layout()
        plt.savefig(f'./plots/curve_validations/{name}.tiff', dpi = 600)
        plt.close()

    #final fitstats dataframe
    df_all.to_csv('curveValidation_withfits.csv')

'''compare splicing outcomes between RNAseq and RT_PCR'''
if args.part2 == True:
    #load df
    df_all = pandas.read_csv('curveValidation_withfits.csv')
    #add deltaPSI columns
    df_all['deltaPSI'] = df_all.apply(lambda x: abs(x['top'] - x['bottom']), axis = 1)
    #adds indexers to merge two dataframes
    df_all['indexer'] = df_all.apply(lambda x: (f'{x.geneSymbol}:{x.exonStart_0base}-{x.exonEnd}'), axis = 1)
   
    #split dataframe on method
    df_1 = df_all[df_all['method'] == 'RNAseq'].reset_index().add_suffix('_RNAseq')
    df_2 = df_all[df_all['method'] == 'RT-PCR']
    #remove unnecessary columns
    df_2 = df_2.drop(['GeneID', 'geneSymbol', 'chr', 'strand', 'exonStart_0base', 'exonEnd',
       'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE', 'exon_size', 'expected_incl_size', 'expected_exc_size'], axis = 1)
    df_2 = df_2.reset_index().add_suffix('_RT-PCR')

    #concat two frames
    df_all = pandas.concat([df_1, df_2], axis = 1)

    ##get read counts
     #get readcount data since I hard coded the RT stuff
    df_orig = pandas.read_pickle('MBNLeventlist_original_HEK.pkl')
    df_orig['indexer_RNAseq'] = df_orig.apply(lambda x: f'{x.geneSymbol}:{x.exonStart_0base}-{x.exonEnd}', axis = 1)
    #only want readcounts
    df_orig = df_orig[['indexer_RNAseq','IJC_SAMPLE_1_d0', 'SJC_SAMPLE_1_d0', 'IJC_SAMPLE_2_d100', 'SJC_SAMPLE_2_d100', 'IJC_SAMPLE_2_d170', 'SJC_SAMPLE_2_d170',
                        'IJC_SAMPLE_2_d250', 'SJC_SAMPLE_2_d250', 'IJC_SAMPLE_2_d280', 'SJC_SAMPLE_2_d280', 'IJC_SAMPLE_2_d1000', 'SJC_SAMPLE_2_d1000',
                        'IJC_SAMPLE_1_avg_d0', 'SJC_SAMPLE_1_avg_d0', 'IJC_SAMPLE_2_avg_d100', 'SJC_SAMPLE_2_avg_d100',
                        'IJC_SAMPLE_2_avg_d170', 'SJC_SAMPLE_2_avg_d170', 'IJC_SAMPLE_2_avg_d250', 'SJC_SAMPLE_2_avg_d250',
                        'IJC_SAMPLE_2_avg_d280', 'SJC_SAMPLE_2_avg_d280', 'IJC_SAMPLE_2_avg_d1000', 'SJC_SAMPLE_2_avg_d1000']]

    #merge dfs
    df_all = df_all.merge(df_orig, how = 'inner', on = 'indexer_RNAseq')
    #custom cols to drop out of laziness
    df_all = df_all.drop(['Unnamed: 0_RNAseq', 'method_RNAseq', 'indexer_RNAseq', 'index_RNAseq', 'index_RT-PCR',
                           'Unnamed: 0_RT-PCR', 'method_RT-PCR', 'indexer_RT-PCR'], axis = 1)
    
   
    #needs wonky conversion to get int from the lists which stored as a string from file to run in this script
    def convert_to_list(s):
        if isinstance(s, str):
            s = s[1:-2]
            s = s.split(',')
            s = [float(x) for x in s]
            return s
    
    #conversion
    cols = ['IncLevel_d0_RNAseq', 'IncLevel_d100_RNAseq', 'IncLevel_d170_RNAseq', 'IncLevel_d250_RNAseq', 'IncLevel_d280_RNAseq', 'IncLevel_d1000_RNAseq',
            'IncLevel_d0_RT-PCR', 'IncLevel_d100_RT-PCR', 'IncLevel_d170_RT-PCR', 'IncLevel_d250_RT-PCR', 'IncLevel_d280_RT-PCR', 'IncLevel_d1000_RT-PCR']
    for col in cols:
        df_all[col] = df_all.apply(lambda x: convert_to_list(x[col]), axis = 1)

    '''plot correlations, all events included'''
    if args.a == True:
        ####################################################################################
        #                   plot slope and LogEC50 correlation
        ####################################################################################
        ####################################### plot slope correlations #################################
        #find the min/max value
        min, max = df_all.loc[:, ['slope_RNAseq', 'slope_RT-PCR']].min().min(), df_all.loc[:, ['slope_RNAseq', 'slope_RT-PCR']].max().max()
        #plot data
        sns.scatterplot(data = df_all, x = 'slope_RNAseq', y = 'slope_RT-PCR')
        #obtain m (slope) and b(intercept) of linear regression line
        m, b, r_value, p_value, std_err = scipy.stats.linregress(df_all['slope_RNAseq'], df_all['slope_RT-PCR'])
        r_sq = r_value**2
        r_value = round(r_value, 3)
        #pval significance
        if p_value < 0.001:
            p_value = "{:.2e}".format(p_value)
        else:
            p_value = "{:.3f}".format(p_value)
        #add linear regression line to scatterplot
        plt.plot(df_all['slope_RNAseq'], (m * df_all['slope_RNAseq'] + b), color = 'black')
        #annotate reg line and stuff
        m, b, r_sq = "{:.2f}".format(m), "{:.2f}".format(b), "{:.2f}".format(r_sq)
        fit_annot = f'y = {m}x+{b}\nR = {r_value}\np-value = {p_value}'
        bbox = dict(edgecolor = 'black', facecolor = 'grey', alpha = 0.5, boxstyle = 'square,pad=0.25') #sets variable to keep annotate flags clean
        plt.annotate(fit_annot, xy = (1.5, 11.5), xycoords = 'data', bbox = bbox, fontsize = 13, horizontalalignment = 'left', verticalalignment = 'center')
        plt.xlabel('Slope RNAseq', fontsize = 'xx-large')
        plt.ylabel('Slope RT-PCR', fontsize = 'xx-large')
        #set limits for graph
        plt.ylim(min - 0.5, max + 0.5)
        plt.xlim(min - 0.5, max + 0.5)
        plt.xticks(fontsize = 17)
        plt.yticks(fontsize = 17)
        plt.title('Slope', fontsize = 'xx-large')
        plt.tight_layout()
        plt.savefig(f'./plots/curve_validations/correlations/RNAseq_vs_RT_slope.tiff', dpi = 600)
        plt.close()
        
        
        ############################### plot LogEC50 correlation #######################################
        #find the min/max value
        min, max = df_all.loc[:, ['LogEC50_RNAseq', 'LogEC50_RT-PCR']].min().min(), df_all.loc[:, ['LogEC50_RNAseq', 'LogEC50_RT-PCR']].max().max()
        #plot data 
        sns.scatterplot(data = df_all, x = 'LogEC50_RNAseq', y = 'LogEC50_RT-PCR')
        #obtain m (slope) and b(intercept) of linear regression line
        m, b, r_value, p_value, std_err = scipy.stats.linregress(df_all['LogEC50_RNAseq'], df_all['LogEC50_RT-PCR'])
        r_sq = r_value**2
        r_value = round(r_value, 3)
        #pval significance
        if p_value < 0.001:
            p_value = "{:.2e}".format(p_value)
        else:
            p_value = "{:.3f}".format(p_value)
        #add linear regression line to scatterplot
        plt.plot(df_all['LogEC50_RNAseq'], (m * df_all['LogEC50_RNAseq'] + b), color = 'black')
        #annotate reg line and stuff
        m, b, r_sq = "{:.2f}".format(m), "{:.2f}".format(b), "{:.2f}".format(r_sq)
        fit_annot = f'y = {m}x+{b}\nR = {r_value}\np-value = {p_value}'
        bbox = dict(edgecolor = 'black', facecolor = 'grey', alpha = 0.5, boxstyle = 'square,pad=0.25') #sets variable to keep annotate flags clean
        plt.annotate(fit_annot, xy = (0.15, 0.84), xycoords = 'data', bbox = bbox, fontsize = 13, horizontalalignment = 'left', verticalalignment = 'center')
        plt.xlabel('LogEC50 RNAseq', fontsize = 'xx-large')
        plt.ylabel('LogEC50 RT-PCR', fontsize = 'xx-large')
        #set limits for graph
        plt.ylim(min - 0.05, max + 0.05)
        plt.xlim(min - 0.05, max + 0.05)
        plt.xticks(fontsize = 17)
        plt.yticks(fontsize = 17)
        plt.title('LogEC50', fontsize = 'xx-large')
        plt.tight_layout()
        plt.savefig(f'./plots/curve_validations/correlations/RNAseq_vs_RT_LogEC50.tiff', dpi = 600)
        plt.close()

        ######################################################################################
        #####                                  plots PSI correlations
        ######################################################################################

        def correlation_scatter(id_vars, hue, cbar_title, labels_only, name, plot_r):
            #reorganize dataframe to make xy plot possible
            val_cols_seq = ['IncLevel_avg_d0_RNAseq', 'IncLevel_avg_d100_RNAseq', 'IncLevel_avg_d170_RNAseq', 'IncLevel_avg_d250_RNAseq', 'IncLevel_avg_d280_RNAseq', 'IncLevel_avg_d1000_RNAseq']
            val_cols_rt = ['IncLevel_avg_d0_RT-PCR', 'IncLevel_avg_d100_RT-PCR', 'IncLevel_avg_d170_RT-PCR', 'IncLevel_avg_d250_RT-PCR', 'IncLevel_avg_d280_RT-PCR', 'IncLevel_avg_d1000_RT-PCR']
            df1 = pandas.melt(df_all, id_vars = id_vars, value_vars = val_cols_seq, var_name = 'variable_RNAseq', value_name = 'value_RNAseq')
            df2 = pandas.melt(df_all, id_vars = id_vars, value_vars = val_cols_rt, var_name = 'variable_RT-PCR', value_name = 'value_RT-PCR')
            df_xy = pandas.concat([df1, df2], axis = 1) #merges dfs
            df_xy = df_xy.loc[:,~df_xy.columns.duplicated()] #drops duplicate deltaPSI column
        
            #this if statement plots the event labels instead of points, know which points belong to which event
            if labels_only == True:
                #label datapoints
                for i in range(len(df_xy)):
                    # plt.text(x = df_xy['value_RNAseq'][i], y = df_xy['value_RT-PCR'][i], s = df_xy['geneSymbol_RNAseq'][i],
                    #         fontdict = dict(color = 'black', size = 5))
                    plt.text(x = df_xy['value_RNAseq'][i], y = df_xy['value_RT-PCR'][i], s = df_xy['geneSymbol_RNAseq'][i],
                            fontdict = dict(color = 'black', size = 5),
                            bbox = dict(boxstyle = 'round4', pad = 0.3, facecolor = 'white', alpha = 1.0))
        
            #this else is the actual plots with points not labels
            else:
                #if statement to void the issue of None when making a colorbar
                if hue != None:
                    norm = plt.Normalize(df_xy[hue].min(), df_xy[hue].max())
                    #can set "Blues_r" to reverse color, make sure to adjust in scatter below
                    sm = plt.cm.ScalarMappable(cmap = "Blues", norm = norm)
                    sm.set_array([])
                    sns.scatterplot(data = df_xy, x = 'value_RNAseq', y = 'value_RT-PCR',
                                     hue = hue, palette = "Blues", edgecolor = 'gray', linewidth = 0.25)
                    plt.colorbar(sm).set_label(cbar_title, fontsize = 15)
                    plt.legend().remove()
        
                #in the case of None, plot
                else:
                    #plot itself
                    sns.scatterplot(data = df_xy, x = 'value_RNAseq', y = 'value_RT-PCR', hue = hue)
        
                #if statement to include regression fit
                if plot_r == True:
                    #obtain m (slope) and b(intercept) of linear regression line
                    m, b, r_value, p_value, std_err = scipy.stats.linregress(df_xy['value_RNAseq'], df_xy['value_RT-PCR'])
                    r_sq = r_value**2
                    r_value = round(r_value, 3)
                    #pval significance
                    if p_value < 0.001:
                        p_value = "{:.2e}".format(p_value)
                    else:
                        p_value = "{:.3f}".format(p_value)
                    #add linear regression line to scatterplot
                    plt.plot(df_xy['value_RNAseq'], (m * df_xy['value_RNAseq'] + b), color = 'black')
                    #annotate reg line and stuff
                    m, b, r_sq = "{:.2f}".format(m), "{:.2f}".format(b), "{:.2f}".format(r_sq) #formats to two decimal places
                    fit_annot = f'y = {m}x+{b}\nR = {r_value}\np-value = {p_value}'
                    bbox = dict(edgecolor = 'black', facecolor = 'grey', alpha = 0.5, boxstyle = 'square,pad=0.25') #sets variable to keep annotate flags clean
                    plt.annotate(fit_annot, xy = (2, 89), xycoords = 'data', bbox = bbox, fontsize = 13, horizontalalignment = 'left', verticalalignment = 'center')
        
            plt.ylim(0, 100)
            plt.xlim(0, 100)
            plt.xlabel('PSI RNAseq', fontsize = 'xx-large')
            plt.ylabel('PSI RT-PCR', fontsize = 'xx-large')
            plt.xticks(fontsize = 17)
            plt.yticks(fontsize = 17)
            plt.tight_layout()
            plt.savefig(f'./plots/curve_validations/correlations/{name}.tiff', dpi = 600)
            plt.close()
        
        correlation_scatter(['deltaPSI_RNAseq', 'geneSymbol_RNAseq'], None, None, False, 'RNAseq_vs_RT_PSIall', True)
        correlation_scatter(['deltaPSI_RNAseq', 'geneSymbol_RNAseq'], 'deltaPSI_RNAseq', None, True, 'RNAseq_vs_RT_PSIall_labelsOnly', True)
        correlation_scatter(['exon_size_RNAseq', 'geneSymbol_RNAseq'], 'exon_size_RNAseq', 'Exon Size', False, 'RNAseq_vs_RT_PSIall_exonHue', True)


        ########################################################################
        ###### correlate exon size and PSI diff ###########
        #########################################################################
        id_vars = ['exon_size_RNAseq', 'geneSymbol_RNAseq']
        #reorganize dataframe to make xy plot possible
        val_cols_seq = ['IncLevel_avg_d0_RNAseq', 'IncLevel_avg_d100_RNAseq', 'IncLevel_avg_d170_RNAseq', 'IncLevel_avg_d250_RNAseq', 'IncLevel_avg_d280_RNAseq', 'IncLevel_avg_d1000_RNAseq']
        val_cols_rt = ['IncLevel_avg_d0_RT-PCR', 'IncLevel_avg_d100_RT-PCR', 'IncLevel_avg_d170_RT-PCR', 'IncLevel_avg_d250_RT-PCR', 'IncLevel_avg_d280_RT-PCR', 'IncLevel_avg_d1000_RT-PCR']
        df1 = pandas.melt(df_all, id_vars = id_vars, value_vars = val_cols_seq, var_name = 'variable_RNAseq', value_name = 'value_RNAseq')
        df2 = pandas.melt(df_all, id_vars = id_vars, value_vars = val_cols_rt, var_name = 'variable_RT-PCR', value_name = 'value_RT-PCR')
        df_xy = pandas.concat([df1, df2], axis = 1) #merges dfs
        df_xy = df_xy.loc[:,~df_xy.columns.duplicated()] #drops duplicate deltaPSI column
        df_xy['PSI_diff_abs'] = df_xy.apply(lambda x: (x['value_RNAseq'] - x['value_RT-PCR']), axis = 1)

        #plot call
        sns.scatterplot(data = df_xy, x = 'exon_size_RNAseq', y = 'PSI_diff_abs')
        #regression line
        #obtain m (slope) and b(intercept) of linear regression line
        m, b, r_value, p_value, std_err = scipy.stats.linregress(df_xy['exon_size_RNAseq'], df_xy['PSI_diff_abs'])
        r_sq = r_value**2
        r_value = round(r_value, 3)
        #pval significance
        if p_value < 0.001:
            p_value = "{:.2e}".format(p_value)
        else:
            p_value = "{:.3f}".format(p_value)
        #add linear regression line to scatterplot
        plt.plot(df_xy['exon_size_RNAseq'], (m * df_xy['exon_size_RNAseq'] + b), color = 'black')
        #annotate reg line and stuff
        m, b, r_sq = "{:.2f}".format(m), "{:.2f}".format(b), "{:.2f}".format(r_sq) #formats to two decimal places
        fit_annot = f'y = {m}x+{b}\nR = {r_value}\np-value = {p_value}'
        bbox = dict(edgecolor = 'black', facecolor = 'grey', alpha = 0.5, boxstyle = 'square,pad=0.25') #sets variable to keep annotate flags clean
        plt.annotate(fit_annot, xy = (8, 42), xycoords = 'data', bbox = bbox, fontsize = 13, horizontalalignment = 'left', verticalalignment = 'center')

        #final changes
        plt.ylim(-30, 50)
        plt.xlim(0, 400)
        plt.xlabel('Exon Size', fontsize = 'xx-large')
        plt.ylabel('PSI RNAseq - PSI RT-PCR', fontsize = 'xx-large')
        plt.xticks(fontsize = 17)
        plt.yticks(fontsize = 17)
        plt.tight_layout()
        plt.savefig(f'./plots/curve_validations/correlations/yExonSize_vs_xPSIdiff.tiff', dpi = 600)
        plt.close()

        ########################################################################
        ###### correlate reads and PSI/PSI diff ###########
        #########################################################################

        df1 = df_all[['IncLevel_d0_RNAseq', 'IncLevel_d0_RT-PCR', 'IncLevel_avg_d0_RNAseq', 'IncLevel_avg_d0_RT-PCR', 'IJC_SAMPLE_1_d0', 'SJC_SAMPLE_1_d0']]
        df2 = df_all[['IncLevel_d100_RNAseq', 'IncLevel_d100_RT-PCR', 'IncLevel_avg_d100_RNAseq', 'IncLevel_avg_d100_RT-PCR', 'IJC_SAMPLE_2_d100', 'SJC_SAMPLE_2_d100']]
        df3 = df_all[['IncLevel_d170_RNAseq', 'IncLevel_d170_RT-PCR', 'IncLevel_avg_d170_RNAseq', 'IncLevel_avg_d170_RT-PCR', 'IJC_SAMPLE_2_d170', 'SJC_SAMPLE_2_d170']]
        df4 = df_all[['IncLevel_d250_RNAseq', 'IncLevel_d250_RT-PCR', 'IncLevel_avg_d250_RNAseq', 'IncLevel_avg_d250_RT-PCR', 'IJC_SAMPLE_2_d250', 'SJC_SAMPLE_2_d250']]
        df5 = df_all[['IncLevel_d280_RNAseq', 'IncLevel_d280_RT-PCR', 'IncLevel_avg_d280_RNAseq', 'IncLevel_avg_d280_RT-PCR', 'IJC_SAMPLE_2_d280', 'SJC_SAMPLE_2_d280']]
        df6 = df_all[['IncLevel_d1000_RNAseq', 'IncLevel_d100_RT-PCR', 'IncLevel_avg_d100_RNAseq', 'IncLevel_avg_d100_RT-PCR', 'IJC_SAMPLE_2_d100', 'SJC_SAMPLE_2_d1000']]

        def adjust_name(df, val):
            df = df.rename(columns={df.columns[0]:'PSI_RNAseq', df.columns[1]:'PSI_RT-PCR', df.columns[2]:'PSIavg_RNAseq', df.columns[3]:'PSIavg_RT-PCR', df.columns[4]:'IJC', df.columns[5]:'SJC'})
            df['dose'] = [[val, val, val]] * len(df)
            return df

        df1 = adjust_name(df1, 'd0')
        df2 = adjust_name(df2, 'd100')
        df3 = adjust_name(df3, 'd170')
        df4 = adjust_name(df4, 'd250')
        df5 = adjust_name(df5, 'd280')
        df6 = adjust_name(df6, 'd1000')
        
        #concats dataframe
        df_reads = pandas.concat([df1, df2, df3, df4, df5, df6])
        
        
        '''
        this section makes plots where the IJC+SJC is averaged and comapred to avg PSI
        '''

        #function with sums IJC and SJC from each replicate, then averages those sums for total read avg
        def avg_reads(l1, l2):
            s1 = l1[0] + l2[0]
            s2 = l1[1] + l2[1]
            s3 = l1[2] + l2[2]
            return (s1 + s2 + s3)/3
        
        df_reads['ijc+sjc_avg'] = df_reads.apply(lambda x: avg_reads(x['IJC'], x['SJC']), axis = 1)
        df_reads['PSIavg_seq-rt'] = df_reads['PSIavg_RNAseq'] - df_reads['PSIavg_RT-PCR']
        df_reads['PSIavg_abs_seq-rt'] = numpy.abs(df_reads['PSIavg_RNAseq'] - df_reads['PSIavg_RT-PCR'])
       
         ########### PSI correlations with reads hue ###############################
        #colorbar and plot
        norm = plt.Normalize(df_reads['ijc+sjc_avg'].min(), df_reads['ijc+sjc_avg'].max())
        #can set "Blues_r" to reverse color, make sure to adjust in scatter below
        sm = plt.cm.ScalarMappable(cmap = "Blues", norm = norm)
        sm.set_array([])
        sns.scatterplot(data = df_reads, x = 'PSIavg_RNAseq', y = 'PSIavg_RT-PCR',
                            hue = 'ijc+sjc_avg', palette = "Blues", edgecolor = 'gray', linewidth = 0.25)
        plt.colorbar(sm).set_label('IJC+SJC Average', fontsize = 15)
        plt.legend().remove()
        #obtain m (slope) and b(intercept) of linear regression line
        m, b, r_value, p_value, std_err = scipy.stats.linregress(df_reads['PSIavg_RNAseq'], df_reads['PSIavg_RT-PCR'])
        r_sq = r_value**2
        r_value = round(r_value, 3)
        #pval significance
        if p_value < 0.001:
            p_value = "{:.2e}".format(p_value)
        else:
            p_value = "{:.3f}".format(p_value)
        #add linear regression line to scatterplot
        plt.plot(df_reads['PSIavg_RNAseq'], (m * df_reads['PSIavg_RNAseq'] + b), color = 'black')
        #annotate reg line and stuff
        m, b, r_sq = "{:.2f}".format(m), "{:.2f}".format(b), "{:.2f}".format(r_sq) #formats to two decimal places
        fit_annot = f'y = {m}x+{b}\nR = {r_value}\np-value = {p_value}'
        bbox = dict(edgecolor = 'black', facecolor = 'grey', alpha = 0.5, boxstyle = 'square,pad=0.25') #sets variable to keep annotate flags clean
        plt.annotate(fit_annot, xy = (2, 89), xycoords = 'data', bbox = bbox, fontsize = 13, horizontalalignment = 'left', verticalalignment = 'center')
        #final adjustments
        plt.ylim(0, 100)
        plt.xlim(0, 100)
        plt.xlabel('PSI RNAseq', fontsize = 'xx-large')
        plt.ylabel('PSI RT-PCR', fontsize = 'xx-large')
        plt.xticks(fontsize = 17)
        plt.yticks(fontsize = 17)
        plt.tight_layout()
        plt.savefig(f'./plots/curve_validations/correlations/RNAseq_vs_RT_PSIall_avg_ReadsHue.tiff', dpi = 600)
        plt.close()

        ############### reads vs PSI abs ###########
        #####plot call for left plot
        sns.scatterplot(data = df_reads, x = 'ijc+sjc_avg', y = 'PSIavg_abs_seq-rt')
        plt.xscale('log')
        plt.xlabel('IJC+SJC Average', fontsize = 'xx-large')
        plt.ylabel('|PSI RNAseq - PSI RT-PCR|', fontsize = 'xx-large')
        plt.xticks(fontsize = 17)
        plt.yticks(fontsize = 17)
        plt.tight_layout()
        plt.savefig(f'./plots/curve_validations/correlations/yreads_vs_xPSIdiffABSs_avg.tiff', dpi = 600)
        plt.close()
        
'''bioinformatic analyses of RT stuff'''
if args.part3 == True:
    '''plot the general distribution of curve params and swarmplots for inc/exc'''
    if args.a == True:
        df_rt = pandas.read_pickle('df_rt_final_pickle.pkl')

        ###### slope and ec50 dist all ########
        #vertical orentation
        #generate a dummy column
        df_rt['dummy'] = ''
        #generate subplot grid
        fig, ax = plt.subplots(1, 2, constrained_layout=True)
        sns.swarmplot(ax = ax[0], x = 'dummy', y = 'slope', data = df_rt, color = 'gray', size = 8)
        sns.boxplot(ax = ax[0], data=df_rt, x='dummy', y='slope', showmeans=True, meanline=True, meanprops={'color': 'k', 'ls': '-', 'lw': 2},
            medianprops={'visible': False}, whiskerprops={'visible': False}, zorder=10,
            showfliers=False, showbox=False, showcaps=False)
        sns.swarmplot(ax = ax[1], x = 'dummy', y = 'LogEC50', data = df_rt, color = 'gray', size = 8)
        sns.boxplot(ax = ax[1], data=df_rt, x='dummy', y='LogEC50', showmeans=True, meanline=True, meanprops={'color': 'k', 'ls': '-', 'lw': 2},
            medianprops={'visible': False}, whiskerprops={'visible': False}, zorder=10,
            showfliers=False, showbox=False, showcaps=False)
        ax[0].set_xlabel('')
        ax[1].set_xlabel('')
        ax[0].set_ylabel('RT-PCR Slope', fontsize = 15)
        ax[1].set_ylabel('RT-PCR LogEC50', fontsize = 15)
        ax[0].tick_params(axis = 'both', which = 'major', labelsize = 15)
        ax[1].tick_params(axis = 'both', which = 'major', labelsize = 15)
        plt.savefig('./plots/RT_only/curvefit_rtOnly_swarmAll_vertical.tiff', dpi = 600)
        plt.close()

        #horizontal orientation
        #generate subplot grid
        fig, ax = plt.subplots(2, 1, constrained_layout=True)
        sns.swarmplot(ax = ax[0],  x = 'slope', y = 'dummy', data = df_rt, color = 'gray', size = 8, orient = 'h')
        sns.boxplot(ax = ax[0], data=df_rt, x ='slope', y ='dummy',  orient = 'h', showmeans=True, meanline=True, meanprops={'color': 'k', 'ls': '-', 'lw': 2},
            medianprops={'visible': False}, whiskerprops={'visible': False}, zorder=10,
            showfliers=False, showbox=False, showcaps=False)
        sns.swarmplot(ax = ax[1], x = 'LogEC50', y = 'dummy',  data = df_rt, color = 'gray', size = 8, orient = 'h')
        sns.boxplot(ax = ax[1], data=df_rt, x ='LogEC50', y ='dummy',  orient = 'h', showmeans=True, meanline=True, meanprops={'color': 'k', 'ls': '-', 'lw': 2},
            medianprops={'visible': False}, whiskerprops={'visible': False}, zorder=10,
            showfliers=False, showbox=False, showcaps=False)
        ax[0].set_xlabel('')
        ax[1].set_xlabel('')
        ax[0].set_ylabel('RT-PCR Slope', fontsize = 15)
        ax[1].set_ylabel('RT-PCR LogEC50', fontsize = 15)
        ax[0].tick_params(axis = 'both', which = 'major', labelsize = 15)
        ax[1].tick_params(axis = 'both', which = 'major', labelsize = 15)
        plt.savefig('./plots/RT_only/curvefit_rtOnly_swarmAll_horizontal.tiff', dpi = 600)
        plt.close()

        '''compares inclusion and exclusion events via swarmplot, due to sample size assume non-parametric tests'''
        def compareIncExc(y, y_ax_title, title):
            #grab values for inclusion and exclusion products
            IncList = (df_rt[y][df_rt['IncLevelDifference'] == 'Inclusion']).reset_index(drop = True).rename('Inclusion')
            ExcList = (df_rt[y][df_rt['IncLevelDifference'] == 'Exclusion']).reset_index(drop = True).rename('Exclusion')

            ####non-parametric test
            test = 'Mann Whitney U Test'
            #tests
            out = scipy.stats.mannwhitneyu(x = IncList, y = ExcList, alternative = 'two-sided')
            pval_2side = out[1]
            out = scipy.stats.mannwhitneyu(x = IncList, y = ExcList, alternative = 'greater')
            pval_greater = out[1]
            out = scipy.stats.mannwhitneyu(x = IncList, y = ExcList, alternative = 'less')
            pval_less = out[1]
            #finding the best p value
            pval_dict = {'2-sided':pval_2side, 'Greater':pval_greater, 'Less':pval_less}
            lowest_p_key = min(pval_dict, key = pval_dict.get)

            ##### plot swarm plot
            #initial plot with some formatting
            df = pandas.concat([IncList, ExcList], axis = 1)
            sns.swarmplot(data = df, size = 10)
            plt.ylabel(y_ax_title, fontsize = 22)
            #set title
            plt.title(f'{test} \n p-value({lowest_p_key}): {pval_dict[lowest_p_key]:.6f}', fontsize = 'small', pad = 20)
            #plot median line
            sns.boxplot(showmeans=False,
                        meanline=False,
                        medianprops={'visible': True, 'color': 'k', 'ls': '-', 'lw': 3},
                        whiskerprops={'visible': False}, zorder=10,
                        data=df, showfliers=False, showbox=False, showcaps=False)
            plt.xticks(fontsize = 22)
            plt.yticks(fontsize = 22)
            plt.rc('axes', labelsize = 22)
            plt.tight_layout()
            plt.savefig('./plots/RT_only/stripplots/non-parametric/%s_swarm.tiff' % title, dpi = 600)
            plt.close()
          
        #removes top and tight spines
        matplotlib.rcParams["axes.spines.right"] = False
        matplotlib.rcParams["axes.spines.top"] = False
        
        #plots
        compareIncExc('slope', 'RT-PCR Slope', 'IncLevel_vs_slope')
        compareIncExc('LogEC50', 'RT-PCR LogEC50', 'IncLevel_vs_LogEC50')
  
'''some RNAseq only stuff'''
if args.part4 == True:
    #this is taken from MBNL_HEK_pipeline_cleaned.py part4, it contains all the RNAseq events
    df_all = pandas.read_pickle('MBNLeventlist_curvefit_HEK.pkl')
    
    '''general comparison between inclusion and exclusion events'''
    if args.a == True:
        #same function with different dataframe from --part4 -b
        def compareIncExc(y, y_ax_title, title):
                #grab values for inclusion and exclusion products
                IncList = (df_all[y][df_all['IncLevelDifference'] == 'Inclusion']).reset_index(drop = True).rename('Inclusion')
                ExcList = (df_all[y][df_all['IncLevelDifference'] == 'Exclusion']).reset_index(drop = True).rename('Exclusion')

                ####non-parametric test
                test = 'Mann Whitney U Test'
                #tests
                out = scipy.stats.mannwhitneyu(x = IncList, y = ExcList, alternative = 'two-sided')
                pval_2side = out[1]
                out = scipy.stats.mannwhitneyu(x = IncList, y = ExcList, alternative = 'greater')
                pval_greater = out[1]
                out = scipy.stats.mannwhitneyu(x = IncList, y = ExcList, alternative = 'less')
                pval_less = out[1]
                #finding the best p value
                pval_dict = {'2-sided':pval_2side, 'Greater':pval_greater, 'Less':pval_less}
                lowest_p_key = min(pval_dict, key = pval_dict.get)

                ##### plot swarm plot
                #initial plot with some formatting
                df = pandas.concat([IncList, ExcList], axis = 1)
                sns.swarmplot(data = df, size = 5)
                plt.ylabel(y_ax_title, fontsize = 22)
                #set title
                plt.title(f'{test} \n p-value({lowest_p_key}): {pval_dict[lowest_p_key]:.6f}', fontsize = 'small', pad = 20)
                #plot median line
                sns.boxplot(showmeans=False,
                            meanline=False,
                            medianprops={'visible': True, 'color': 'k', 'ls': '-', 'lw': 3},
                            whiskerprops={'visible': False}, zorder=10,
                            data=df, showfliers=False, showbox=False, showcaps=False)
                plt.xticks(fontsize = 22)
                plt.yticks(fontsize = 22)
                plt.rc('axes', labelsize = 22)
                plt.tight_layout()
                plt.savefig('./plots/swarmplots/%s_swarm.tiff' % title, dpi = 600)
                plt.close()

        #removes top and tight spines
        matplotlib.rcParams["axes.spines.right"] = False
        matplotlib.rcParams["axes.spines.top"] = False 

        #plots
        compareIncExc('LogEC50', 'RNAseq LogEC50', 'IncLevel_vs_LogEC50_stat')
        compareIncExc('slope', 'RNAseq Slope', 'IncLevel_vs_slope_stat')

    
