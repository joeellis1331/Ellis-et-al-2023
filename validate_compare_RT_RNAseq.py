#!/usr/bin/python3

import argparse
import numpy
import pandas
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import math
import scipy
from scipy.optimize import curve_fit
from Bio import SeqIO
from collections import Counter
from maxentpy import maxent
from maxentpy.maxent import load_matrix5, load_matrix3

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
'''
if args.part2 == True:

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
        xdata = [math.log10(1), math.log10(2), math.log10(2.7), math.log10(3.5), math.log10(3.8), math.log10(11)]
        ydata = [row.IncLevel_avg_d0, row.IncLevel_avg_d100, row.IncLevel_avg_d170, row.IncLevel_avg_d250, row.IncLevel_avg_d280, row.IncLevel_avg_d1000]

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

        #define plus minus symbol for plotting
        plus_minus = (u'\u00B1')

        #RNAseq data
        if df_all.iloc[ind].method == 'RNAseq':
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
            plt.plot(x, y,  c = 'orange', label = '%s fit' % df_all.iloc[ind].method) #plots sigmoid fit

            #appends parameters for table plotting to array
            top_table = f'{round(max([popt[0], popt[1]]), 1)}{plus_minus}{df_all.iloc[ind].stderr_top}'
            bottom_table = f'{round(min([popt[0], popt[1]]), 1)}{plus_minus}{df_all.iloc[ind].stderr_bottom}'
            ec50_table = f'{round(popt[3], 2)}{plus_minus}{df_all.iloc[ind].stderr_LogEC50}'
            slope_table = f'{round(popt[2], 2)}{plus_minus}{df_all.iloc[ind].stderr_slope}'
            r_sq_table = f'{round(df_all.iloc[ind].R_sq, 3)}'
            extend_list = [top_table, bottom_table, ec50_table, slope_table, r_sq_table]
            cellText[0].extend(extend_list)

        #RT-PCR data
        else:
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
            plt.plot(x, y,  c = 'purple', label = '%s fit' % df_all.iloc[ind].method) #plots sigmoid fit

            #appends parameters for table plotting to array
            top_table = f'{round(max([popt[0], popt[1]]), 1)}{plus_minus}{df_all.iloc[ind].stderr_top}'
            bottom_table = f'{round(min([popt[0], popt[1]]), 1)}{plus_minus}{df_all.iloc[ind].stderr_bottom}'
            ec50_table = f'{round(popt[3], 2)}{plus_minus}{df_all.iloc[ind].stderr_LogEC50}'
            slope_table = f'{round(popt[2], 2)}{plus_minus}{df_all.iloc[ind].stderr_slope}'
            r_sq_table = f'{round(df_all.iloc[ind].R_sq, 3)}'
            extend_list = [top_table, bottom_table, ec50_table, slope_table, r_sq_table]
            cellText[1].extend(extend_list)

        name = '%s_%s-%s' % (df_all.iloc[ind].geneSymbol, str(df_all.iloc[ind].exonStart_0base), str(df_all.iloc[ind].exonEnd))
        return name, cellText


    #necessary plotting stuff, goes by tuple index
    ind = [x for x in range(len(df_all))] #creates list of indexes for len of dataframe
    it = iter(ind) #creates iterator
    ind = list(zip(it, it)) #accesses iterator and advances position, creating a tuple of pairs e.g. ((1,2), (3,4))
   
    #iterates through index tuple
    for tup in ind:
        #2D list for plotting table
        cellText = [['RNAseq'], ['RT-PCR']]
        #plots and grabs geneSymbol for chart title
        name, cellText = plot_curves(tup[0], cellText)
        name, cellText = plot_curves(tup[1], cellText)

        #plots table
        colLabels = ['$\\bf{Method}$', '$\\bf{Top}$', '$\\bf{Bottom}$', '$\\bf{Slope}$', '$\\bf{LogEC50}$', '$\\bf{R\u00b2}$'] #uses TeX code
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
        plt.title(df_all.iloc[tup[0]].geneSymbol, pad = 54, fontsize = 'xx-large', style = 'italic')
        plt.tight_layout()
        plt.savefig(f'./plots/curve_validations/{name}.tiff', dpi = 600, format = 'png')
        plt.close()

    #final fitstats dataframe
    df_all.to_csv('curveValidation_withfits_original.csv')

'''compare splicing outcomes between RNAseq and RT_PCR'''
if args.part3 == True:
    df_all = pandas.read_csv('curveValidation_withfits_original.csv')

    #needs wonky conversion to get int from the lists which stored as a string from file to run in this script
    def convert_to_list(s):
        s = s[1:-2]
        s = s.split(',')
        s = [float(x) for x in s]
        return s

    cols = ['IncLevel_d0', 'IncLevel_d100', 'IncLevel_d170', 'IncLevel_d250', 'IncLevel_d280', 'IncLevel_d1000']
    for col in cols:
        df_all[col] = df_all.apply(lambda x: convert_to_list(x[col]), axis = 1)

    #add deltaPSI columns
    df_all['deltaPSI'] = df_all.apply(lambda x: abs(x['top'] - x['bottom']), axis = 1)
    #add no log EC50
    df_all['EC50'] = df_all.apply(lambda x: (10**x['LogEC50']), axis = 1)


    #more data for readcounts
    df_fitstats = pandas.read_pickle('MBNLeventlist_fitstats_HEK.pkl')
    #adds indexers to merge two dataframes
    df_fitstats['indexer'] = df_fitstats.apply(lambda x: f'{x.geneSymbol}:{x.exonStart_0base}-{x.exonEnd}', axis = 1)
    df_all['indexer'] = df_all.apply(lambda x: (f'{x.geneSymbol}:{x.exonStart_0base}-{x.exonEnd}'), axis = 1)
    #only want readcounts
    df_fitstats = df_fitstats[['indexer','IJC_SAMPLE_1_avg_d0', 'SJC_SAMPLE_1_avg_d0', 'IJC_SAMPLE_2_avg_d100', 'SJC_SAMPLE_2_avg_d100',
                                'IJC_SAMPLE_2_avg_d170', 'SJC_SAMPLE_2_avg_d170', 'IJC_SAMPLE_2_avg_d250', 'SJC_SAMPLE_2_avg_d250',
                                'IJC_SAMPLE_2_avg_d280', 'SJC_SAMPLE_2_avg_d280', 'IJC_SAMPLE_2_avg_d1000', 'SJC_SAMPLE_2_avg_d1000']]
    #merge dfs
    df_all = df_all.merge(df_fitstats, how = 'inner', on = 'indexer').drop(columns = 'Unnamed: 0')

    #calculate sum of the avg reads for all samples
    df_all['counts_d0'] = df_all.apply(lambda x: (x.IJC_SAMPLE_1_avg_d0 + x.SJC_SAMPLE_1_avg_d0), axis = 1)
    df_all['counts_d100'] = df_all.apply(lambda x: (x.IJC_SAMPLE_2_avg_d100 + x.SJC_SAMPLE_2_avg_d100), axis = 1)
    df_all['counts_d170'] = df_all.apply(lambda x: (x.IJC_SAMPLE_2_avg_d170 + x.SJC_SAMPLE_2_avg_d170), axis = 1)
    df_all['counts_d250'] = df_all.apply(lambda x: (x.IJC_SAMPLE_2_avg_d250 + x.SJC_SAMPLE_2_avg_d250), axis = 1)
    df_all['counts_d280'] = df_all.apply(lambda x: (x.IJC_SAMPLE_2_avg_d100 + x.SJC_SAMPLE_2_avg_d280), axis = 1)
    df_all['counts_d1000'] = df_all.apply(lambda x: (x.IJC_SAMPLE_2_avg_d1000 + x.SJC_SAMPLE_2_avg_d1000), axis = 1)
    df_all['total_avg_reads'] = df_all.apply(lambda x: (x.IJC_SAMPLE_1_avg_d0 + x.SJC_SAMPLE_1_avg_d0 + x.IJC_SAMPLE_2_avg_d100 + x.SJC_SAMPLE_2_avg_d100 + x.IJC_SAMPLE_2_avg_d170 + x.SJC_SAMPLE_2_avg_d170 + x.IJC_SAMPLE_2_avg_d250 + x.SJC_SAMPLE_2_avg_d250 + x.IJC_SAMPLE_2_avg_d280 + x.SJC_SAMPLE_2_avg_d280 + x.IJC_SAMPLE_2_avg_d1000), axis = 1)
    df_all['log(total_reads)'] = df_all.apply(lambda x: math.log10(x.total_avg_reads), axis = 1)

    #concat the datasets so each even is on a single row
    df_1 = df_all[df_all['method'] == 'RNAseq'].reset_index().add_suffix('_RNAseq')
    df_2 = df_all[df_all['method'] == 'RT-PCR']
    #remove extra columns
    df_2 = df_2[['GeneID', 'geneSymbol', 'chr', 'strand', 'exonStart_0base', 'exonEnd',
       'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE',
       'IncLevel_d0', 'IncLevel_d100', 'IncLevel_d170', 'IncLevel_d250',
       'IncLevel_d280', 'IncLevel_d1000', 'IncLevel_avg_d0',
       'IncLevel_avg_d100', 'IncLevel_avg_d170', 'IncLevel_avg_d250',
       'IncLevel_avg_d280', 'IncLevel_avg_d1000', 'method', 'RSS', 'LogRSS',
       'top', 'bottom', 'LogEC50', 'slope', 'pooled_StDev', 'deltaPSI', 'EC50',
       'indexer']]
    df_2 = df_2.reset_index().add_suffix('_RT-PCR')
    df_all = pandas.concat([df_1, df_2], axis = 1).drop(columns = 'index_RNAseq')

    '''plot correlations, all events included'''
    if args.a == True:
        ####################################################################################
        #                   plot slope and LogEC50 correlation
        ####################################################################################
        ####################################### plot slope correlations #################################
        #find the min/max value
        min, max = df_all.loc[:, ['slope_RNAseq', 'slope_RT-PCR']].min().min(), df_all.loc[:, ['slope_RNAseq', 'slope_RT-PCR']].max().max()
        
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
        plt.annotate(fit_annot, xy = (2.5, 11.5), xycoords = 'data', bbox = bbox, fontsize = 13, horizontalalignment = 'left', verticalalignment = 'center')
        plt.xlabel('Slope RNAseq', fontsize = 'xx-large')
        plt.ylabel('Slope RT-PCR', fontsize = 'xx-large')
        #set limits for graph
        plt.ylim(min + 0.5, max + 0.5)
        plt.xlim(min + 0.5, max + 0.5)
        plt.xticks(fontsize = 17)
        plt.yticks(fontsize = 17)
        plt.title('Slope', fontsize = 'xx-large')
        plt.tight_layout()
        plt.savefig(f'./plots/curve_validations/RNAseq_vs_RT_slope.tiff', dpi = 600, format = 'png')
        plt.close()
        
        
        ############################### plot LogEC50 correlation #######################################
        #find the min/max value
        min, max = df_all.loc[:, ['LogEC50_RNAseq', 'LogEC50_RT-PCR']].min().min(), df_all.loc[:, ['LogEC50_RNAseq', 'LogEC50_RT-PCR']].max().max()
        
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
        plt.annotate(fit_annot, xy = (0.23, 0.89), xycoords = 'data', bbox = bbox, fontsize = 13, horizontalalignment = 'left', verticalalignment = 'center')
        plt.xlabel('LogEC50 RNAseq', fontsize = 'xx-large')
        plt.ylabel('LogEC50 RT-PCR', fontsize = 'xx-large')
        #set limits for graph
        plt.ylim(min + 0.05, max + 0.05)
        plt.xlim(min + 0.05, max + 0.05)
        plt.xticks(fontsize = 17)
        plt.yticks(fontsize = 17)
        plt.title('LogEC50', fontsize = 'xx-large')
        plt.tight_layout()
        plt.savefig(f'./plots/curve_validations/RNAseq_vs_RT_LogEC50.tiff', dpi = 600, format = 'png')
        plt.close()

    
        ######################################################################################
        #####                                  plots PSI correlations
        ######################################################################################

        def correlation_scatter(id_vars, hue, labels_only, name, plot_r):
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
                    sm = plt.cm.ScalarMappable(cmap = "Blues_r", norm = norm)
                    sm.set_array([])
                    sns.scatterplot(data = df_xy, x = 'value_RNAseq', y = 'value_RT-PCR', hue = hue, palette = "Blues_r")
                    title = name.split('_')[-1]
                    plt.colorbar(sm).set_label(title, fontsize = 15)
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
            plt.savefig(f'./plots/curve_validations/{name}.tiff', dpi = 600, format = 'png')
            plt.close()
        
        correlation_scatter(['deltaPSI_RNAseq', 'geneSymbol_RNAseq'], None, False, 'RNAseq_vs_RT_PSIall', True)
        correlation_scatter(['deltaPSI_RNAseq', 'geneSymbol_RNAseq'], 'deltaPSI_RNAseq', True, 'RNAseq_vs_RT_PSIall_labelsOnly', True)

'''bioinformatic analyses of RT stuff'''
if args.part4 == True:
    '''getting the dataframe and all the pre stats ready for analysis'''
    if args.a == True:
        df_all = pandas.read_csv('curveValidation_withfits_original.csv')

        #needs wonky conversion to get int from the lists which stored as a string from file to run in this script
        def convert_to_list(s):
            s = s[1:-2]
            s = s.split(',')
            s = [float(x) for x in s]
            return s

        cols = ['IncLevel_d0', 'IncLevel_d100', 'IncLevel_d170', 'IncLevel_d250', 'IncLevel_d280', 'IncLevel_d1000']
        for col in cols:
            df_all[col] = df_all.apply(lambda x: convert_to_list(x[col]), axis = 1)

        #add deltaPSI columns
        df_all['deltaPSI'] = df_all.apply(lambda x: abs(x['top'] - x['bottom']), axis = 1)

        #more data for readcounts
        df_fitstats = pandas.read_pickle('MBNLeventlist_fitstats_HEK.pkl')
        #adds indexers to merge two dataframes
        df_fitstats['indexer'] = df_fitstats.apply(lambda x: f'{x.geneSymbol}:{x.exonStart_0base}-{x.exonEnd}', axis = 1)
        df_all['indexer'] = df_all.apply(lambda x: (f'{x.geneSymbol}:{x.exonStart_0base}-{x.exonEnd}'), axis = 1)
        #only want readcounts
        df_fitstats = df_fitstats[['indexer','IJC_SAMPLE_1_avg_d0', 'SJC_SAMPLE_1_avg_d0', 'IJC_SAMPLE_2_avg_d100', 'SJC_SAMPLE_2_avg_d100',
                                    'IJC_SAMPLE_2_avg_d170', 'SJC_SAMPLE_2_avg_d170', 'IJC_SAMPLE_2_avg_d250', 'SJC_SAMPLE_2_avg_d250',
                                    'IJC_SAMPLE_2_avg_d280', 'SJC_SAMPLE_2_avg_d280', 'IJC_SAMPLE_2_avg_d1000', 'SJC_SAMPLE_2_avg_d1000']]
        #merge dfs
        df_all = df_all.merge(df_fitstats, how = 'inner', on = 'indexer').drop(columns = 'Unnamed: 0')

        #calculate sum of the avg reads for all samples
        df_all['counts_d0'] = df_all.apply(lambda x: (x.IJC_SAMPLE_1_avg_d0 + x.SJC_SAMPLE_1_avg_d0), axis = 1)
        df_all['counts_d100'] = df_all.apply(lambda x: (x.IJC_SAMPLE_2_avg_d100 + x.SJC_SAMPLE_2_avg_d100), axis = 1)
        df_all['counts_d170'] = df_all.apply(lambda x: (x.IJC_SAMPLE_2_avg_d170 + x.SJC_SAMPLE_2_avg_d170), axis = 1)
        df_all['counts_d250'] = df_all.apply(lambda x: (x.IJC_SAMPLE_2_avg_d250 + x.SJC_SAMPLE_2_avg_d250), axis = 1)
        df_all['counts_d280'] = df_all.apply(lambda x: (x.IJC_SAMPLE_2_avg_d100 + x.SJC_SAMPLE_2_avg_d280), axis = 1)
        df_all['counts_d1000'] = df_all.apply(lambda x: (x.IJC_SAMPLE_2_avg_d1000 + x.SJC_SAMPLE_2_avg_d1000), axis = 1)
        df_all['total_avg_reads'] = df_all.apply(lambda x: (x.IJC_SAMPLE_1_avg_d0 + x.SJC_SAMPLE_1_avg_d0 + x.IJC_SAMPLE_2_avg_d100 + x.SJC_SAMPLE_2_avg_d100 + x.IJC_SAMPLE_2_avg_d170 + x.SJC_SAMPLE_2_avg_d170 + x.IJC_SAMPLE_2_avg_d250 + x.SJC_SAMPLE_2_avg_d250 + x.IJC_SAMPLE_2_avg_d280 + x.SJC_SAMPLE_2_avg_d280 + x.IJC_SAMPLE_2_avg_d1000), axis = 1)
        df_all['log(total_reads)'] = df_all.apply(lambda x: math.log10(x.total_avg_reads), axis = 1)

        #grab just the RT data
        df_rt = df_all[df_all['method'] == 'RT-PCR']
        #remove extra columns
        df_rt = df_rt[['GeneID', 'geneSymbol', 'chr', 'strand', 'exonStart_0base', 'exonEnd',
           'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE',
           'IncLevel_d0', 'IncLevel_d100', 'IncLevel_d170', 'IncLevel_d250',
           'IncLevel_d280', 'IncLevel_d1000', 'IncLevel_avg_d0',
           'IncLevel_avg_d100', 'IncLevel_avg_d170', 'IncLevel_avg_d250',
           'IncLevel_avg_d280', 'IncLevel_avg_d1000', 'method', 'RSS', 'LogRSS',
           'top', 'bottom', 'LogEC50', 'slope', 'pooled_StDev', 'deltaPSI',
           'indexer']]
        df_rt = df_rt.reset_index()

        ########### get sequences ###########
        #function to convert df to a bed file format
        def makeBed(df):
            #up exon coordinates
            #copy of data and saves to new dataframe
            df_upExon = pandas.DataFrame(df[['chr', 'upstreamES', 'upstreamEE', 'geneSymbol', 'strand']].copy())
            #adds blank columns to make strand column in proper place
            df_upExon.insert(4, 'blank', '1')
            #removes column titles
            df_upExon = df_upExon.set_axis([1, 2, 3, 4, 5, 6], axis = 1)

            #upintron coordinates
            df_upIntron = pandas.DataFrame(df[['chr', 'upstreamEE', 'exonStart_0base', 'geneSymbol', 'strand']].copy())
            df_upIntron.insert(4, 'blank', '1')
            df_upIntron = df_upIntron.set_axis([1, 2, 3, 4, 5, 6], axis = 1)
            #AS exon coordinates
            df_exon = pandas.DataFrame(df[['chr', 'exonStart_0base', 'exonEnd', 'geneSymbol', 'strand']].copy())
            df_exon.insert(4, 'blank', '1')
            df_exon = df_exon.set_axis([1, 2, 3, 4, 5, 6], axis = 1)
            #downintron coordinates
            df_downIntron = pandas.DataFrame(df[['chr', 'exonEnd', 'downstreamES', 'geneSymbol', 'strand']].copy())
            df_downIntron.insert(4, 'blank', '1')
            df_downIntron = df_downIntron.set_axis([1, 2, 3, 4, 5, 6], axis = 1)
            #down exon coords
            df_downExon = pandas.DataFrame(df[['chr', 'downstreamES', 'downstreamEE', 'geneSymbol', 'strand']].copy())
            df_downExon.insert(4, 'blank', '1')
            df_downExon = df_downExon.set_axis([1, 2, 3, 4, 5, 6], axis = 1)

            #concats to one dataframe
            bed_df = pandas.concat([df_upExon, df_upIntron, df_exon, df_downIntron, df_downExon])
            #removes the 'chr' prefix of the chr column for use in getFasta from bedtools
            bed_df[1] = bed_df[1].str[3:]
            return bed_df
        #makes bed files
        # bed_out = makeBed(df_all)
        # bed_out.to_csv('./BED_RTevents_HEK.txt', sep = '\t', index = False, header = False)


        #matches fasta sequences back to dataframe
        #appends sequences to dataframe using the fasta dictionary
        def append_seq(df, seq_dict):
            #creates indexes
            upExon_index = '%s:%s-%s' % (df['chr'], df['upstreamES'], df['upstreamEE'])
            upIntron_index = '%s:%s-%s' % (df['chr'], df['upstreamEE'], df['exonStart_0base'])
            exon_index = '%s:%s-%s' % (df['chr'], df['exonStart_0base'], df['exonEnd'])
            downIntron_index = '%s:%s-%s' % (df['chr'], df['exonEnd'], df['downstreamES'])
            downExon_index = '%s:%s-%s' % (df['chr'], df['downstreamES'], df['downstreamEE'])
            #matches sequences based on indexes
            df['upExon'] = seq_dict[upExon_index]
            df['upIntron'] = seq_dict[upIntron_index]
            df['exon'] = seq_dict[exon_index]
            df['downIntron'] = seq_dict[downIntron_index]
            df['downExon'] = seq_dict[downExon_index]
            #extra needed stuff
            df['upIntron250'] = df['upIntron'][-251:-1]
            df['downIntron250'] = df['downIntron'][:250]
            df['upIntron_length'] = len(df['upIntron'])
            df['exon_length'] = len(df['exon'])
            df['downIntron_length'] = len(df['downIntron'])
            df['upIntron_length250'] = len(df['upIntron250'])
            df['downIntron_length250'] = len(df['downIntron250'])
            return df
        #pasrse fasta, wraps seqSplit
        def fasta_parse(input_fasta, df):
            seq_dict = {} #empty dictionary
            with open(input_fasta) as handle:
                #iterate through each fasta record
                for record in SeqIO.FastaIO.FastaIterator(handle): #note: this utilizes a SeqRecord class
                    indexer = 'chr' + str(record.id) #creates match indexer with chr, up and down coords
                    indexer = indexer[:-3] #removes strand parameter
                    seq = str(record.seq) #grabs sequence
                    if indexer not in seq_dict:
                        seq_dict[indexer] = seq #appends to a dictionary
                    else:
                        continue

            #append sequences to dataframe
            df = df.apply(lambda x: append_seq(x, seq_dict), axis = 1)
            #add extra column stuff
            return df
        #parses sequences from fasta file into a dictionary and then matches to dataframe, should be unique
        df_rt = fasta_parse('/mnt/c/Users/joeel/Desktop/Transient_Coreg_Files/NEW_bioinfo/HEK/fasta_RTevents_HEK.txt', df_rt)

        ####### calculate splice site scores using maxentscan py wrapper from https://github.com/kepbod/maxentpy #########
        #function to get SS scores for up and down 5SS and 3SS
        def ss_scores(df):
            #preload matrices
            matrix5 = load_matrix5()
            matrix3 = load_matrix3()
            #get sequences sections
            seq5_up = df['upExon'][-3:] + df['upIntron'][:6]
            seq3_up = df['upIntron'][-20:] + df['exon'][:3]
            seq5_down = df['exon'][-3:] + df['downIntron'][:6]
            seq3_down = df['downIntron'][-20:] + df['downExon'][:3]
            #store the sequence chunks for reference
            df['up_5SS_seq'] = seq5_up
            df['down_5SS_seq'] = seq5_down
            df['up_3SS_seq'] = seq3_up
            df['down_3SS_seq'] = seq3_down
            #run SS scores
            df['up_5SS'] = maxent.score5(seq5_up, matrix = matrix5)
            df['down_5SS'] = maxent.score5(seq5_down, matrix = matrix5)
            df['up_3SS'] = maxent.score3(seq3_up, matrix = matrix3)
            df['down_3SS'] = maxent.score3(seq3_down, matrix = matrix3)
            return df
        #get splice site scores
        df_rt = df_rt.apply(lambda x: ss_scores(x), axis = 1)

        #single value of all splice site scores using methodology from the pub PMC3950716
        #S = (up_5ss - down_5ss) + (down_3SS - up_3SS)
        #higher S-values reflect a higher priority for pairing upI_ss5 and dnI_ss3, leading to the skipping isoform
        df_rt['spliceCompS'] = (df_rt['up_5SS'] - df_rt['down_5SS']) + (df_rt['down_3SS'] - df_rt['up_3SS'])

        ###### finding YGCYs #############
        #suppresses copy warning
        pandas.options.mode.chained_assignment = None  # default='warn'
        '''function taken from part 3 and modified'''
        #function to parse a string into kmers
        def get_k_mer(x, k):
           length = len(x)
           return [x[i: i+ k].lower() for i in range(length-k+1)] #keeps all sequence information lower case

        #function to get motif counts
        def YGCY_counts(x):
            #get counts for kmers
            indv_count = Counter(x)
            #gets actual YGCY total
            y_count = indv_count['tgct'] +  indv_count['cgcc'] + indv_count['tgcc'] + indv_count['cgct']
            return y_count

        #get kmers in 250nt segment of intronic space
        df_rt['kmer_up'] = df_rt.apply(lambda x: get_k_mer(x['upIntron'], 4), axis = 1)
        df_rt['kmer_up250'] = df_rt.apply(lambda x: get_k_mer(x['upIntron250'], 4), axis = 1)
        df_rt['kmer_exon'] = df_rt.apply(lambda x: get_k_mer(x['exon'], 4), axis = 1)
        df_rt['kmer_down'] = df_rt.apply(lambda x: get_k_mer(x['downIntron'], 4), axis = 1)
        df_rt['kmer_down250'] = df_rt.apply(lambda x: get_k_mer(x['downIntron250'], 4), axis = 1)
        #get counts
        df_rt['YGCY_upIntron'] = df_rt.apply(lambda x: YGCY_counts(x['kmer_up']), axis = 1)
        df_rt['YGCY_upIntron250'] = df_rt.apply(lambda x: YGCY_counts(x['kmer_up250']), axis = 1)
        df_rt['YGCY_exon'] = df_rt.apply(lambda x: YGCY_counts(x['kmer_exon']), axis = 1)
        df_rt['YGCY_downIntron'] = df_rt.apply(lambda x: YGCY_counts(x['kmer_down']), axis = 1)
        df_rt['YGCY_downIntron250'] = df_rt.apply(lambda x: YGCY_counts(x['kmer_down250']), axis = 1)
        #get frequency
        df_rt['YGCY_freq_upIntron250'] = df_rt['YGCY_upIntron250']/df_rt['upIntron_length250']
        df_rt['YGCY_freq_exon250'] = df_rt['YGCY_exon']/df_rt['exon_length']
        df_rt['YGCY_freq_downIntron250'] = df_rt['YGCY_downIntron250']/df_rt['downIntron_length250']
        df_rt['YGCY_freq_upInt250_exon'] = (df_rt['YGCY_upIntron250'] + df_rt['YGCY_exon'])/(df_rt['upIntron_length250'] + df_rt['exon_length'])
        df_rt['YGCY_freq_all_250'] = (df_rt['YGCY_upIntron250'] + df_rt['YGCY_exon'] + df_rt['YGCY_downIntron250'])/(df_rt['upIntron_length250'] + df_rt['exon_length'] + df_rt['downIntron_length250'])

        #find YGCY motif positions
        def findMotifs(kmer_list, flip):
            pos_list = []
            for i, kmer in enumerate(kmer_list):
                if kmer == 'tgct':
                    pos_list.append(i)
                elif kmer == 'cgcc':
                    pos_list.append(i)
                elif kmer == 'tgcc':
                    pos_list.append(i)
                elif kmer == 'cgct':
                    pos_list.append(i)
            
            #transposes upintron to negative values
            if flip != True:
                return pos_list
            else:
                return [-(250-x) for x in pos_list]
        
        df_rt['upInt250_YGCYpos'] = df_rt.apply(lambda x: findMotifs(x['kmer_up250'], True), axis = 1)
        df_rt['exon_YGCYpos'] = df_rt.apply(lambda x: findMotifs(x['kmer_exon'], False), axis = 1)
        df_rt['downInt250_YGCYpos'] = df_rt.apply(lambda x: findMotifs(x['kmer_down250'], False), axis = 1)
        
        ##### getting distance of YGCY to exon-intron boundaries ######
        def getDistance(seq, kmer_list):
            YGCY_list = ['tgct', 'cgcc', 'tgcc', 'cgct']
            pos_list = []
            for i, k in enumerate(kmer_list):
                if k in YGCY_list:
                    pos_list.append(i)
            
            if 'up' in seq:
                return [-(250-x) for x in pos_list]
            else:
                return pos_list
        
        df_rt['YGCY_proxy_upIntron250'] = df_rt.apply(lambda x: getDistance('upIntron', x['kmer_up250']), axis = 1)
        df_rt['YGCY_proxy_exon'] = df_rt.apply(lambda x: getDistance('exon', x['kmer_exon']), axis = 1)
        df_rt['YGCY_proxy_downIntron250'] = df_rt.apply(lambda x: getDistance('downIntron', x['kmer_down250']), axis = 1)


        ######## some extra stuff ################
        df_rt['combined_seq'] = df_rt['upIntron'] + df_rt['exon'] + df_rt['downIntron']
        #seq lengths
        df_rt['seq_length_up'] = df_rt.apply(lambda x: len(x['upIntron']), axis = 1 )
        df_rt['seq_length_exon'] = df_rt.apply(lambda x: len(x['exon']), axis = 1 )
        df_rt['seq_length_down'] = df_rt.apply(lambda x: len(x['downIntron']), axis = 1 )
        df_rt['seq_length_all'] = df_rt.apply(lambda x: len(x['combined_seq']), axis = 1 )
        #GC content
        df_rt['GC_content_up'] = df_rt.apply(lambda x:
                    (((x['upIntron'].count('G') + x['upIntron'].count('C')) / x['seq_length_up']) * 100), axis = 1)
        df_rt['GC_content_exon'] = df_rt.apply(lambda x:
                    (((x['exon'].count('G') + x['exon'].count('C')) / x['seq_length_exon']) * 100), axis = 1)
        df_rt['GC_content_down'] = df_rt.apply(lambda x:
                    (((x['downIntron'].count('G') + x['downIntron'].count('C')) / x['seq_length_down']) * 100), axis = 1)
        df_rt['GC_content'] = df_rt.apply(lambda x:
                    (((x['combined_seq'].count('G') + x['combined_seq'].count('C')) / x['seq_length_all']) * 100), axis = 1)

        #this checks the GC content outside of the presence of YGCYs
        def GC_noY_all(df_rt):
            total_YGCY = df_rt['YGCY_upIntron'] + df_rt['YGCY_exon'] + df_rt['YGCY_downIntron']
            total_YGCY = total_YGCY * 2 #since each YGCY has a G and a C
            GC = df_rt['combined_seq'].count('G') + df_rt['combined_seq'].count('C')
            GCnoY = ((GC - total_YGCY) / df_rt['seq_length_all']) * 100
            return GCnoY

        df_rt['GC_content_noY'] = df_rt.apply(lambda x: GC_noY_all(x), axis = 1)

        def GC_noY(df_rt, seqY, seq, seqL):
            total_YGCY = df_rt[seqY]
            total_YGCY = total_YGCY * 2 #since each YGCY has a G and a C
            GC = df_rt[seq].count('G') + df_rt[seq].count('C')
            GCnoY = ((GC - total_YGCY) / df_rt[seqL]) * 100
            return GCnoY

        df_rt['GC_content_noY_up'] = df_rt.apply(lambda x: GC_noY(x, 'YGCY_upIntron', 'upIntron', 'seq_length_up'), axis = 1)
        df_rt['GC_content_noY_exon'] = df_rt.apply(lambda x: GC_noY(x, 'YGCY_exon', 'exon', 'seq_length_exon'), axis = 1)
        df_rt['GC_content_noY_down'] = df_rt.apply(lambda x: GC_noY(x, 'YGCY_downIntron', 'downIntron', 'seq_length_down'), axis = 1)


        #dataframe with inclusion or exclusion events only
        df_rt['IncLevelDifference'] = df_rt['IncLevel_avg_d0'] - df_rt['IncLevel_avg_d1000']
        df_rt['IncLevelDifference'] = numpy.where((df_rt['IncLevelDifference'] < 0.0), 'Inclusion', 'Exclusion')

        #save upandas.ted dataframe for use later
        df_rt.to_pickle('df_rt_final_pickle.pkl')
        df_rt.to_excel('df_rt_final.xlsx')

    '''plot the general distribution of curve params and swarmplots for inc/exc'''
    if args.b == True:
        df_rt = pandas.read_pickle('df_rt_final_pickle.pkl')

        ### extra comparisons
        df_rt['all_150_YGCY'] = df_rt['YGCY_upIntron250'] + df_rt['YGCY_exon'] + df_rt['YGCY_downIntron250']
        df_rt['up250+exon_YGCY'] = df_rt['YGCY_upIntron250'] + df_rt['YGCY_exon']

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
        plt.savefig('./plots/RT_only/curvefit_rtOnly_swarmAll_vertical.tiff', dpi = 600, format = 'png')
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
        plt.savefig('./plots/RT_only/curvefit_rtOnly_swarmAll_horizontal.tiff', dpi = 600, format = 'png')
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
            plt.savefig('./plots/RT_only/stripplots/non-parametric/%s_swarm.tiff' % title, dpi = 600, format = 'png')
            plt.close()
          
        #removes top and tight spines
        matplotlib.rcParams["axes.spines.right"] = False
        matplotlib.rcParams["axes.spines.top"] = False
        
        #plots
        compareIncExc('slope', 'RT-PCR Slope', 'IncLevel_vs_slope')
        compareIncExc('LogEC50', 'RT-PCR LogEC50', 'IncLevel_vs_LogEC50')
  
'''some RNAseq only stuff'''
if args.part5 == True:
    #this is taken from MBNL_HEK_pipeline_cleaned.py part4, it contains all the RNAseq events
    df_all = pandas.read_pickle('correlation_pickle.pkl')
    
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
                plt.savefig('./plots/swarmplots/%s_swarm.tiff' % title, dpi = 600, format = 'png')
                plt.close()

        #removes top and tight spines
        matplotlib.rcParams["axes.spines.right"] = False
        matplotlib.rcParams["axes.spines.top"] = False 

        #plots
        compareIncExc('LogEC50', 'RNAseq LogEC50', 'IncLevel_vs_LogEC50_stat')
        compareIncExc('slope', 'RNAseq Slope', 'IncLevel_vs_slope_stat')

    