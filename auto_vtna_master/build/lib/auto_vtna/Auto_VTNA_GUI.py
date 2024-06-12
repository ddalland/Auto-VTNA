import PySimpleGUI as sg
import threading
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import logging
logging.getLogger('matplotlib.font_manager').disabled = True
import auto_vtna
import tkinter as tk
from auto_vtna.Normal_VTNA import Normal_VTNA
from auto_vtna.VTNA_functions import VTNA_omissions
import copy
import re
import os
import csv
import tkinter as tk
from tkinter import scrolledtext
from auto_vtna.VTNA_functions import VTNA_tT, on_click, VTNA_new, score_fit,data_scaling_y, origin_line, simple_line
import itertools
from functools import partial
from openpyxl.drawing.image import Image
from io import BytesIO
import matplotlib.ticker as tkr
import warnings
import time
from polyfit import PolynomRegressor, Constraints
from num2words import num2words
from scipy.optimize import curve_fit
# Define a global variable to use to cancel calculations if requested by user. 
global cancel_calculation
cancel_calculation=False
print('Welcome to the Automatic VTNA Calculator! Visit https://github.com/ddalland/Auto-VTNA for more details or check out the Auto-VTNA pre-print: https://chemrxiv.org/engage/chemrxiv/article-details/65fddc2d66c1381729948bb2')
print('Developed by Daniel Dalland with help and input from Dr. Linden Schrecker and Prof. King Kuok (Mimi) Hii at the Imperial College Department of Chemistry.')
print('We are grateful for the creators of the following dependencies:')
print('1. PySimpleGUI: Simplified GUI creation.')
print('2. Pandas: Data manipulation and analysis.')
print('3. NumPy: Numerical computing.')
print('4. Matplotlib: Data visualization.')
print('5. kinter: GUI toolkit.')
print('6. Openpyxl: Excel file handling.')
print('7. Polyfit: Monotonically constrained polynomial fitting.')
print('8. Num2words: Number to word conversion.')
print('10. Scipy: Curve fitting.')
print("We are also grateful for OpenAI's ChatGTP which was used when writing parts of the code.")
print("For more information, check out the Auto-VTNA Github page: https://github.com/ddalland/Auto-VTNA and the pre-print of the relevant publication: https://chemrxiv.org/engage/chemrxiv/article-details/66269b0321291e5d1d64c3c9")

def extract_range_info(values,data):
    '''Function used in the range mode data cropping window to extract 
    inputted value ranges.
    Arguments: 
     - values: data croppting range mode window inputs.
     - data: kinetic data dictionary.'''
    RMS_ss={}
    # Loop through all experiments in the dataset.
    for i in list(data.keys())+['total']:
        start_stop_values=[]
        flag=0
        for j in list(values.keys()):
            # Collect the range values for a given RM.
            if i in j and 'less' not in j and len(values[j])>0 and j.replace(i,'')[0]=='_':
                flag=1
                start_stop_values.append(values[j])
        # If range values were found, add the list to the new dict.
        if flag==1:
            RMS_ss[i]=start_stop_values
    return RMS_ss

def check_ends(values,data):
    ''' Function used in the range mode data cropping window to check
    that range value inputs are correct.
    Arguments: 
     - values: data croppting range mode window inputs.
     - data: kinetic data dictionary.'''
    # Returns correct as True if the end and start values have been 
    # loaded pairwise correctly. 
    vals=[]
    correct=True
    values_2=values.copy()
    # Remove non START and STOP values from the values dict.
    for i in values.keys():
        if 'START' not in i and 'END' not in i:
            values_2.pop(i)
    # Delete the key:value pairs without any value.
    for i in values.keys():
        if len(values[i])==0:
            del values_2[i]
    # Check that start and end values are present pairwise.
    for i in list(data.keys())+['total']:
        value=0
        for j in values_2.keys():
            if i in j:
                value+=1
        vals.append(value)
    # If any 1 or 3 values are present, this means that not all pairs are complete. 
    if 1 in vals or 3 in vals:
        correct=False
        return correct
    # Find the end of the keys in the values_2 dict.
    ends=[i[-1] for i in list(values_2.keys())]
    # Check that the number of ends is an even number.
    if len(ends)%2!=0:
        correct=False
        return correct
    # Check again for pairwise values.
    for i in range(0,len(ends),2):
        if ends[i]!=ends[i+1]:
            correct=False
    return correct

def check_increase(range_info):
    ''' Function used in the range mode data cropping window to check
    that range value pairs are increasing. Avoids errors in the auto_vtna 
    function VTNA_omissions()
    Argument: range_info - output from extract_range_info() '''
    correct=True
    for i in range_info.keys():
        for j in range(0,len(range_info[i]),2):
            if float(range_info[i][j])>=float(range_info[i][j+1]):
                correct=False
    return correct

def is_float(string):
    ''' Exception handling function used to check if a string can be 
    converted to a float without giving an error.'''
    try:
        float_value = float(string)
        return True
    except ValueError:
        return False

def overlay_plot(kinetic_data,experiments,NS_dict,output_reaction_species,tt_scaler=1,y_unit='M',\
                 size_scaler=1,grid=False,title=None,fit_metric='RMSE',constraint='Monotonic',\
                 deg=5,show_fit_function=False,show_overlay_score=False,DP_scaler=1,xtick_rotation=0,\
                 show_legend=True,legend_outside=False,extra_legend=True,line_scaler=1,
                 tT_notation='Normal'):
    ''' Creates a VTNA_selection dictionary for the GUI settings stored in the arguments
    "experiments", "NS_dict" and "output_reaction_species", calculates the transformed time axis 
    for the kinetic data "kinetic_data" and generates the overlay plot for the selected plot settings.
    '''
    # Generate a VTNA selection dictionary template.
    Default_VTNA_selection_dictionary=auto_vtna.make_VTNA_selection(kinetic_data)
    VTNA_selection_dictionary={}
    # Generate a sub-dictionary to specify which experiments to include based on sheet selections.
    VTNA_selection_dictionary['RM']={}
    for i in experiments:
        VTNA_selection_dictionary['RM'][i]=Default_VTNA_selection_dictionary['RM'][i]
    # Update the normalised species and output species for generating the overlay plot based on GUI selections.
    VTNA_selection_dictionary['normalised_species']=NS_dict
    VTNA_selection_dictionary['output_species']=output_reaction_species
    # Calculate the normalised time axis using Normal_VTNA from auto_vtna.
    Run_VTNA=Normal_VTNA(kinetic_data,VTNA_selection_dictionary)
    if constraint=='Via origin':
        constraint='origin'
    elif constraint=='Monotonic':
        constraint='monotonic'
    # Generate the overlay plot.
    Run_VTNA.plot_VTNA_with_overlay_score(y_unit=y_unit,size_scaler=size_scaler,grid=grid, \
                  xaxis_size_scaler=tt_scaler,title=title,fit_metric=fit_metric,constraint=constraint,\
                  deg=deg,show_fit_function=show_fit_function,show_overlay_score=show_overlay_score,\
                    DP_scaler=DP_scaler,xtick_rotation=xtick_rotation,show_legend=show_legend,
                    legend_outside=legend_outside,extra_legend=extra_legend,linewidth=line_scaler,
                    xaxis_notation=tT_notation)
class Automatic_VTNA:
    def __init__(self,data,VTNA_selection,order_range=[-1.5,2.5],resolution=7,deg=5, \
                 fit_metric='RMSE',iterations=7,constraint='monotonic',score_interval=0.15,\
                fixed_order_species=None,initial_mesh_denser=True):
        # Initialise the class:
        self.data=data.copy()
        self.VTNA_selection=VTNA_selection.copy()
        self.order_range=order_range
        self.resolution=resolution
        self.deg=deg
        self.fit_metric=fit_metric
        self.iterations=iterations
        self.constraint=constraint
        self.score_interval=score_interval
        self.fixed_order_species=fixed_order_species
        self.initial_mesh_denser=initial_mesh_denser
        # Perform automatic VTNA calculation:
        Calculation_results = VTNA_orders(self.data,self.VTNA_selection, \
        self.order_range,self.resolution,self.deg,self.fit_metric,self.iterations,self.constraint, \
        self.score_interval,self.fixed_order_species,self.initial_mesh_denser)
        if Calculation_results!=None:
            self.results, self.best_orders, self.interval, self.best_slope=Calculation_results
        else:
            self.results=None

    def identify_error_intervals(self,score_interval=0.15):
        """
        Generates a Pandas dataframe containing the lowest and highest order value which is found in the \
        list of order value combinations that give an overlay score within "score_interval" of the optimal overlay \
        score. The score_interval can be set to any value above 0. 
        """
        # Get the column names except the rightmost one
        col_names = self.results.columns[:-1]
        # Create an empty list to store lists of column values
        combined_values = [[] for _ in range(len(self.results))]
        # Iterate through each row of the DataFrame
        for index, row in self.results.iterrows():
            # Iterate through each column except the rightmost one
            for col_name in col_names:
                # Append the value of each column to the corresponding list
                combined_values[index].append(row[col_name])
        # Convert the list of lists into a numpy 2D array
        order_combinations = np.array(combined_values)
        # Extract the rightmost column as a numpy array
        FQs = self.results.iloc[:, -1].values
        # Obtain the minimum point order combination (minium FQ for variance, RMSE and SE, max FQ for R2).
        minFQ=min(FQs)
        if self.fit_metric=='R2':
            minFQ=max(FQs)
        # Identify order combinations with overlay scores up to "score_interval" higher than that of the best order:
        if self.fit_metric=='R2':
            good_order_combos=[order_combinations[i] for i in range(len(order_combinations)) if FQs[i]>minFQ*(1-score_interval)]
        else: 
            good_order_combos=[order_combinations[i] for i in range(len(order_combinations)) if FQs[i]<minFQ*(1+score_interval)]
        # Identify the smallest and largest order values in the good order combinations list for each normalised species.
        largest_values=[max(column) for column in zip(*good_order_combos)]
        smallest_values=[min(column) for column in zip(*good_order_combos)]
        intervals=list(zip(smallest_values, largest_values))
        intervals=np.array(list(zip(smallest_values, largest_values)))
        # Define a list of the normalised reaction species that haven't been fixed. 
        if type(self.fixed_order_species)!=list:
            fixed_order_list=[self.fixed_order_species]
        else:
            fixed_order_list=self.fixed_order_species
        # Define the list of normalised species and generate the overlay score interval dictionary. 
        normalised_species=[i for i in list(self.VTNA_selection['normalised_species'].keys()) if i not in fixed_order_list]
        interval_dict={'Normalised species':normalised_species,f'Lower order limit {score_interval*100}%':intervals[:,0],f'Upper order limit {score_interval*100}%':intervals[:,1]}
        return pd.DataFrame(interval_dict)

    def plot_orders_vs_overlay(self,points=False,size_scaler=1,size_scaler2=1,title=None,zoom=False,grid=False,color_scaler=1,
    fixed_cbar_resolution=None,interval=False,y_unit='M',title_fontsize=12,overlay_score_range_max=None,
    overlay_score_range_min=None,show_legend_main=True,show_legend_popup=True,decimals_in_colorbar=3,annotate_optimum=True,
    specified_score_interval=None):
        """
        Creates a graph or contour plot visualising the overlay score versus order value matrix \
        for one and two normalised reaction species respectively. If the overlay score versus \
        order value matrix involves 3 or more normalised species, no plot is created. The graph \
        or contour plots can be clicked upon to generate the VTNA overlay plot corresponding to \
        the order values of the click. 
        Args:
            points (bol): Set to True to visualise the order value data points for which overlay \
            scores have been determined. 
            size_scaler (int or float): Factor by which the size of the order value versus overlay \
            score plot is adjusted. 
            size_scaler2 (int or float): Factor by which the size of the VTNA overlay plot generated \
            by clicking on the original plot is adjusted. 
            title (str): The title of the order value versus overlay score plot.
            zoom (int or float): Width of order values around the optimial order point that will be \
            included in the order value versus overlay score plot.
            grid (bol): Set to True to include a grid in the order value versus overlay score plot.
            color_scaler (int or float, ∈<0,1]): Number between 1 and 0 which determines the factor by \
            which the contour plot color bar range is truncated. NB: This will give white sections that \
            are off the contour plot color map.
            resolution_contour (float): Difference in overlay score between separate color shades in \
            the contour plot. 
            interval (bol): Set to True to show each order combination datapoint within the interval_score \
            times the best overlay score value away from the optimal order value point. 
            y_unit (str): The concentration unit of the output concentration profiles. Gets included in the \
            yaxis title. 
            title_fontsize (int): Fontsize of the title of the order value versus overlay score plot.
            overlay_score_range_max (None or float): The maximum overlay score value to be shown in the plot \
            either as a contour plot color bar (2 NS) or the y axis of an ordinary graph (1 NS).
            overlay_score_range_min (None or float): The minimum overlay score value to be shown in the plot \
            either as a contour plot color bar (2 NS) or the y axis of an ordinary graph (1 NS).
            show_legend_main (bol): Set to True to show the legend of the overlay score versus order plot. 
            show_legend_popup (bol): Set to True to show the legend of the VTNA overlay plots generated by \
            mouseclicks on the overlay score versus order plot. 
            decimals_in_colorbar (int): The number of decimals to be shown in the colorbar of an overlay score \
            versus order contour plot. 
            annotate_optimum (bol): Set to True to show the order values of the optimum point of an overlay score \
            versus order contour plot. 

        """
        #if color_scaler>1 or color_scaler<=0:
        #    raise ValueError('color_scaler should be set to ∈<0,1]')
        # Find the species corresponding to the x or x and y axis:
        axis_species=[i.replace('order in ','') for i in list(self.results.columns)[:-1]]
        output_species=self.VTNA_selection['output_species']
        analysis_columns=self.results.columns
        no_normalised_species=len(analysis_columns)-1
        # Create string to include information about reaction species with fixed order value(s) during the 
        # automatic VTNA calculation.
        fixed_order_string=''
        if self.fixed_order_species!=None:
            if type(self.fixed_order_species)==str:
                fixed_order_string=fixed_order_string+f"\nFixed order for {self.fixed_order_species}: {self.VTNA_selection['normalised_species'][self.fixed_order_species]}"
            else:
                for i in self.fixed_order_species:
                    fixed_order_string=fixed_order_string+f"\nFixed order for {i}: {self.VTNA_selection['normalised_species'][i]}."
        if interval==True:
            # Identify the score interval from the calculation or from the specified_score_interval argument.
            score_interval=self.score_interval
            if specified_score_interval!=None:
                if float(specified_score_interval)>0:
                    score_interval=float(specified_score_interval)
                else:
                    warnings.warn("The specified score interval must be a number above 0. The calculation value will be used instead.")
        if no_normalised_species>2:
            return ValueError('The overlay score : order value matrix can only be visualized for 1 or 2 varying normalised reaction species.')
        if no_normalised_species==1:
            x_best=self.best_orders[0]
            # Define x values as orders and y values as fit metric values.
            x_s=np.array(self.results[analysis_columns[0]])
            y_s=np.array(self.results[analysis_columns[1]])
            # Sort the order values from smallest to largest to get one continuous plt line. 
            xy=list(itertools.zip_longest(x_s,y_s))
            xy.sort()
            # Redefine sorted x and y values
            y=[y[1] for y in xy]
            x=[x[0] for x in xy]
            # Find the finest orders interval used to find minimum point
            resolution_orders=str(format(x_s[-1]-x_s[-2],".12f"))
            # Find the index for first non 0 digit in resolution
            for i in range(len(resolution_orders)):
                if resolution_orders[i]!='.' and resolution_orders[i]!='0':
                    break
            # Increase i by 1 if first character in order is '-'
            if str(x_best)[0]=='-':
                i=i+1
            # Create figure with size scaled by size_scaler
            fig, ax = plt.subplots(figsize=(6*size_scaler,4*size_scaler))
            # Plot the fit metric values obtained at each order.
            ax.plot(x,y,marker='o',linestyle='-',markersize=5)
            # Plot the minimum point and insert label information for legend.
            minimum_info=f'Best order: {f"%.{i-1}f" % round(x_best,i-1)}\nUncertainty: {f"%.{i}f" % round(float(resolution_orders[:(i+2)]),i)}\nMinimum {self.fit_metric}: {f"%.{i+2}f" % round(min(y),i+2)}'
            optimum_value=min(y)
            if self.fit_metric=='R2':
                minimum_info=f'Best order: {f"%.{i-1}f" % round(x_best,i-1)}\nUncertainty: {f"%.{i}f" % round(float(resolution_orders[:(i+2)]),i)}\nMaximum {self.fit_metric}: {f"%.{i+2}f" % round(max(y),i+2)}'
                optimum_value=max(y)
            minimum_info=minimum_info+fixed_order_string
            plt.plot(x_best,optimum_value,color='orange',marker='o',label=minimum_info)
            # Define the x and y axis labels. 
            ax.set_xlabel(f'{analysis_columns[0]}',fontsize=12)
            ax.set_ylabel(f'Overlay score: {self.fit_metric}',fontsize=12)
            # Change the x- and y-axis limits in case zoom has been applied. 
            if zoom!=False:
                x1=x_best-zoom/2
                x2=x_best+zoom/2
                plt.xlim(x1,x2)
                y2=[i for i in y if x[y.index(i)]<x2 and x[y.index(i)]>x1]
                y1=min(y2)
                y2=max(y2)
                leeway=(y2-y1)*0.1
                ax.set_ylim(y1-leeway,y2+leeway)
            # Plot a horizontal line showing the datapoints within the interval_score*best overlay score range of the minimum point.
            if interval==True:
                interval_limit=min(y)+min(y)*score_interval
                if self.fit_metric=='R2':
                    interval_limit=max(y)-max(y)*score_interval
                ax.plot([self.order_range[0],self.order_range[1]],[interval_limit,interval_limit],linestyle='--',label=f"{float(score_interval)*100}% from optimal point.")
            # Set the figure title. 
            if title==None:
                fig.suptitle(f'VTNA overlay score versus orders for the formation of {output_species}',fontsize=title_fontsize,y=0.95)
            else:
                fig.suptitle(title,fontsize=title_fontsize,y=0.95)
            # Apply legend and grid if selected.
            if show_legend_main:
                ax.legend()
            if grid==True:
                plt.grid()
            # Define the ylimits if overlay score range max or min have been defined. 
            if overlay_score_range_max!=None:
                ax.set_ylim(top=overlay_score_range_max)
            if overlay_score_range_min!=None:
                ax.set_ylim(bottom=overlay_score_range_min)
        # Code for contour plot to illustrate the overlay score across order values if there are 2 normalised species. 
        if no_normalised_species==2:
            # Obtain the order value arrays x_s and y_s and the overlay scores z_s.
            z_s=np.array(self.results[analysis_columns[2]])
            x_s=np.array(self.results[analysis_columns[0]])
            y_s=np.array(self.results[analysis_columns[1]])
            # Define the x and y values of the minimum point (or maximum for R2).
            x_best=self.best_orders[0]
            y_best=self.best_orders[1]
            # Create figure with size scaled by size_scaler.
            fig, ax = plt.subplots(figsize=(6*size_scaler,5*size_scaler))
            # Identify the max and min fit metric values for setting up color grid for contour plot.
            min_z=min(z_s)
            max_z=max(z_s)
            # Find the finest orders interval used to find minimum point.
            resolution_orders=str(format(y_s[-1]-y_s[-2],".12f"))
            # Find the index for first non 0 digit in resolution_orders.
            for i in range(len(resolution_orders)):
                if resolution_orders[i]!='.' and resolution_orders[i]!='0':
                    break
            # Define q and w to determine number of digits to show for each minimum point order.
            q,w=i,i
            # Increase q and w by 1 if first character in respective order is '-'
            if str(x_best)[0]=='-':
                q=1+i
            if str(y_best)[-1]=='-':
                w=1+i
            if points==True:
                # Plot the data points used to generate the countour plot.
                ax.plot(x_s,y_s,'or',markersize=0.5)
            # Plot the datapoints with an overlay score off the minimum point (or max for R2) score.
            
            if interval==True:
                if self.fit_metric!='R2':
                    x_interval=[x_s[i] for i in range(len(x_s)) if z_s[i]<(score_interval+1)*min(z_s)]
                    y_interval=[y_s[i] for i in range(len(y_s)) if z_s[i]<(score_interval+1)*min(z_s)]
                else:
                    x_interval=[x_s[i] for i in range(len(x_s)) if z_s[i]>(-score_interval+1)*min(z_s)]
                    y_interval=[y_s[i] for i in range(len(y_s)) if z_s[i]>(-score_interval+1)*min(z_s)]
                # Plot the identified interval datapoints in green.
                ax.plot(x_interval,y_interval,markersize=0.7,linestyle='',marker='x',color='b',  
                        label=f'Orders within {score_interval*100}% of optimal O.S.',zorder=2)
            # Print axis labels with species_order from data to assign normalised species.
            ax.set_xlabel(f'{analysis_columns[0]}')
            ax.set_ylabel(f'{analysis_columns[1]}')
            # Plot the data point corresponding to the best VTNA alignment and use label to provide information.
            if self.fit_metric=='R2':
                ax.plot(x_best,y_best,marker='x',color='k',linestyle='',label = f'Best orders: ({round(x_best, q-1):.{int(q-1)}f},{round(y_best, w-1):.{int(w-1)}f})\nUncertainty: {resolution_orders[:(i+2)]}\nMaximum {self.fit_metric}: {round(max(z_s), i+2):.{int(i+2)}f}{fixed_order_string}')
            else:
                ax.plot(x_best,y_best,marker='x',color='k',linestyle='',label=f'Best orders: ({round(x_best, q-1):.{int(q-1)}f},{round(y_best, w-1):.{int(w-1)}f})\nUncertainty: {resolution_orders[:(i+2)]}\nMinimum {self.fit_metric}: {round(min(z_s), i+2):.{int(i+2)}f}{fixed_order_string}')
            # If zoom argument applied, apply xlim and ylim around the minimum point and find the lowest z_value in that interval.
            min_z=min(z_s)
            max_z=max(z_s)
            z_range=max_z-min_z
            max_z=max(z_s)-z_range*(1-color_scaler)
            # If the overlay score range max or min arguments have been defined, prepare the relevant variables.
            cbar_fixing=False
            if overlay_score_range_max!=None:
                cbar_fixing=True
                max_z=overlay_score_range_max
            if overlay_score_range_min!=None:
                cbar_fixing=True
                min_z=overlay_score_range_min
            if zoom!=False:
                if cbar_fixing==True:
                    warnings.warn("overlay_score_range_max and overlay_score_range_min is ignored when zoom isn't set to False.")
                # Define x and y limits to zoom in closer to the minimum point
                min_z=min(z_s)
                ax.set_xlim(x_best-zoom/2,x_best+zoom/2)
                ax.set_ylim(y_best-zoom/2,y_best+zoom/2)
                # Obtain the z_s values in the new domain.
                z_s2=[i for i in z_s if x_s[list(z_s).index(i)]<(x_best+zoom/2) and x_s[list(z_s).index(i)]>(x_best-zoom/2) \
                       and y_s[list(z_s).index(i)]<(y_best+zoom/2) and y_s[list(z_s).index(i)]>(y_best-zoom/2)] 
                # Define the maximum fit metric value in this interval and the difference between this and the minimum value
                max_z=max(z_s2)
                z_range=max_z-min_z
                # Use color_scaler to define to adjust the maximum value of the colorbar
                max_z=max_z-z_range*(1-color_scaler)
            # Define the order value range of the plot
            delta_order=max(x_s)-min(x_s)
            # Define appropriate colorbar values if R2 is used as fit metric
            if self.fit_metric=='R2':
                max_z=max(z_s)
                max_z=max_z+(1-max_z)*0.2
                min_z=min_z+z_range*(1-color_scaler)
            if fixed_cbar_resolution!=None:
                resolution_cbar=float(fixed_cbar_resolution)
            else:
                resolution_cbar=abs((max_z-min_z)/30)
            clev = np.arange(min_z,max_z,resolution_cbar)
            # Create contour plot.
            CS = ax.tricontourf(x_s, y_s, z_s, clev,cmap=plt.cm.viridis_r)
            if self.fit_metric=='R2':
                CS = ax.tricontourf(x_s, y_s, z_s, clev,cmap=plt.cm.viridis)
            cbar=fig.colorbar(CS,orientation='vertical',label=f'Overlay score: {self.fit_metric}',ax=ax,format=tkr.FormatStrFormatter(f'%.{decimals_in_colorbar}g'))
            legend_loc='best'
            if zoom==False:
                best_value_x=x_best+0.045
                best_value_y=y_best+0.045
                if best_value_y-min(y_s)<best_value_y-max(y_s):
                    legend_loc='upper '
                else: 
                    legend_loc='lower '
                if best_value_x-min(x_s)<best_value_x-max(x_s):
                    legend_loc=legend_loc+'right'
                else:
                    legend_loc=legend_loc+'left'
            if show_legend_main:
                ax.legend(loc=legend_loc)
            # Define the title of the figure.
            if title==None:
                fig.suptitle(f'VTNA overlay score versus orders for the formation of {output_species}',fontsize=title_fontsize,y=0.95)
            else:
                fig.suptitle(title,fontsize=title_fontsize,y=0.95)
            # Annotate the minimum point to show its order values. 
            # Adjust position according to the "zoom" value and delta_order if zoom isn't set to False
            if annotate_optimum:
                if zoom!=False:
                    ax.annotate(f'{f"%.{q}f" % round(x_best,q)}, {f"%.{w}f" % round(y_best,w)}',(x_best+0.045*zoom/delta_order, y_best+0.045*zoom/delta_order))
                if zoom==False:
                    ax.annotate(f'{f"%.{q}f" % round(x_best,q)}, {f"%.{w}f" % round(y_best,w)}',(best_value_x,best_value_y))
        # Create a new click function with inputs.
        on_click_with_inputs = partial(on_click,data=self.data,VTNA_selection=self.VTNA_selection,y_unit=y_unit,\
        size_scaler=size_scaler2,fit_metric=self.fit_metric,deg=self.deg,constraint=self.constraint,fixed_order_species=self.fixed_order_species,
        axis_species=axis_species,show_legend=show_legend_popup)
        # Bind the on_click function to the mouse click event.
        cid = fig.canvas.mpl_connect('button_press_event', on_click_with_inputs)
        # Show the plot.
        plt.tight_layout()
        if grid==True:
            plt.grid()
        image_stream = BytesIO()
        plt.savefig(image_stream, format='png', bbox_inches='tight', pad_inches=0)
        image_stream.seek(0)
        img = Image(image_stream)
        plt.show()
        return img

def VTNA_orders(data,VTNA_selection,order_range=[-1.5,2.5],resolution=8,deg=5,fit_metric='RMSE',iterations=6,constraint='monotonic',
                score_interval=0.15,fixed_order_species=None,initial_mesh_denser=True):
        """
        Calculates the overlay score across different reaction order values for each normalised species in the \
        VTNA_selection dictionary for the kinetic data "data". First, a uniform order grid is defined by the order_range and \
        resolution arguments. Then, VTNA_omissions is applied to remove experiments or datapoints within experiments \
        as defined in VTNA_selection. Next, the average output species concentration profiles are normalised to \
        1. Then, numpy arrays containing the data required for variable time normalisation using the trapezoid \
        rule is defined with padding to ensure compatibility with numpy and inputed to explore_orders() to calculate \
        the overlay score for each order combination by calculating the normalised time axis with VTNA_tT() \
        and performing a global fit of each time normalised output concentration profile with a fitting method defined by \
        keyword arguments "deg" and "constraint". Explore_orders() then returns the overlay scores, each order combination \
        investigated and the best order combination overall. Next, another order grid with values around each minimum order value \
        is defined and explored by explore_orders() which calculates the overlay score for each combination. The new minimum \
        point is then used to define an even more fine order combination grid which is again given to explore_orders(). This \
        is repeated "iterations" number of time to improve the precision of the assigned order values. In the end, a \
        pandas DataFrame conatining the order values for each normalised species and the resulting overlay scores is \
        returned with the best order values and the intervals within which the overlay score deviates from the optimum \
        value by a factor of less than the "score_interval" value and the best order values found.
        Args:
            data (dict): Kinetic data as a dictionary of Pandas dataframes.
            VTNA_selection (dict): Dictionary containing information about how to carry out automatic VTNA.
            order_range (list): Minimum and maximum order value used to define the initial order value grid. 
            resolution (int): Number of datapoints included in the order value grids in each explore_orders() iteration.
            deg (int): Polynomial degree used in the global fitting prodecure to calculate the overlay score. 
            fit_metric (str): Goodness-of-fit metric chosen to define the overlay score value. Either 'R2', 'RMSE' or 'SE'. \
            otherwise variance is used.
            iterations (int): The number of times explore_orders() is called after the first order grid has been investigated.
            constraint (str): Defines the constraint of the global fit used by explore_orders() to calculate the overlay score. \
            Can be set to 'monotonic' for monotonic fitting (default), None for ordinary polynomial fitting or 'origin' for \
            linear fit through origin.
            score_interval: Defines the cutoff minimum_overlay_score + minimum_overlay_score*score_interval (max for R2) for \
            quoting the range of order values that give overlay scores close to that of the best order values to quantify how \
            well defined best order values are. 
            fixed_order_species (list or None): List containing the normalised reaction species for which \
            time axis should be normalised with a fixed value (as defined in the VTNA selection dictionary). \
            This can be done to lower the dimensionality of the overlay score versus order matrix. 
            initial_mesh_denser (bol): Determines whether the initial order value mesh should be denser than the \
            resolution argument by a factor of 1.5 (True) or the same as the resolution parameter (False).
        """
        # Check that the constraint argument is specified correctly.
        if constraint!=None and constraint!='origin' and constraint!='monotonic':
            raise ValueError("Constraint can only be set to 'origin' for linear fit through origin, 'monotonic' for \
                          monotonic fit or None for ordinary polynomial fit")
        # Check the intial time:
        start_time=time.time()
        print('Auto-VTNA calculation underway.')
        # Check if the number of fixed order species is greater than the number of normalised reaction species.
        if fixed_order_species!=None:
            no_fixed_order_species=len(fixed_order_species)
            if type(fixed_order_species)==str:
                no_fixed_order_species=1
            if no_fixed_order_species>=len(VTNA_selection['normalised_species'].items()):
                raise ValueError("The number of fixed order species can't be the same as or lower than the \
    number of normalised reaction species in the VTNA selection dictionary.")
        # Obtain deepcopy of the kinetic data and VTNA selection dictionary to avoid alterations. 
        data_copy=copy.deepcopy(data)
        selection_copy=copy.deepcopy(VTNA_selection)
        # Access one of the experimental dataframes.
        example_RM=list(data_copy.keys())[0]
        # Define the time axis key as t.
        t=data_copy[example_RM].keys()[0] 
        # Use the selected resolution and order range to define a list of order combinations. 
        # Resultion is increased by 50% for the first order exploration.
        if initial_mesh_denser:
            resolution_scaler=1.5
        else:
            resolution_scaler=1
        orders=np.linspace(order_range[0],order_range[1],round(resolution*resolution_scaler))
        # Edit the kinetic data according to the VTNA selection dictionary and add a scaled output values column.
        data_omissions=VTNA_omissions(data_copy,selection_copy)
        data_scaled=data_scaling_y(data_omissions,selection_copy)
        # Collect the time, normalised species and output species concentration data into designated lists.
        output_lists=[]
        NS_lists=[]
        time_lists=[]
        normalised_species=list(selection_copy['normalised_species'].keys())
        # Check that the reaction species listed in 'fixed_order_species' are present as normalised reaction 
        # species in the VTNA selection dictionary. Remove these reaction species from the normalised_species list. 
        if fixed_order_species!=None:
            # Make fixed_order_species into a list if only one reaction species string is defined. 
            if type(fixed_order_species)==str:
                fixed_order_species=[fixed_order_species]
            # Cancel calculation if fixed order reaction species have been selected that are not defined
            # as normalised reaction species in the VTNA selection dictionary. 
            if all(species in normalised_species for species in fixed_order_species)!=True:
                raise ValueError("All reaction species listed in 'fixed_order_species' must be present as a normalised \
                species in the selection dictionary")
            # Remove the fixed reaction species from the normalised_species list.
            normalised_species=[i for i in normalised_species if i not in fixed_order_species]
            # Generate a list of the fixed order values from the VTNA selection dictionary. 
            fixed_orders=[selection_copy['normalised_species'][i] for i in fixed_order_species]
            NS_lists_fixed_orders=[]
        # Find the number of normalised species and add a duplicate of orders to order_combinations for each.
        no_orders=len(normalised_species)
        order_combinations=[]
        for i in range(no_orders):
            order_combinations.append(orders)
        # Obtain the names of the experiment datasets in the data dictionary.
        experiments=list(data_copy.keys())
        # If selected experiments are defined in the selection dictionary, use these experiment names instead.
        if 'RM' in list(selection_copy.keys()):
            experiments=list(selection_copy['RM'].keys())
        # Extract the information required for time variable normalisation and fitting across order combinations. 
        counter=0
        for i in normalised_species:
            NS_nested_list=[]
            for exp in experiments:
                # Obtain the concentration profiles for each of the normalised reaction species i from each experiment exp.
                NS_nested_list.append(list(data_scaled[str(exp)][i].to_numpy()))
                # Collect the time and scaled output concentration values for each experiment once. 
                if counter==0:
                    time_lists.append(list(data_scaled[str(exp)][t].to_numpy()))
                    output_lists.append(list(data_scaled[str(exp)]['output_scaled'].to_numpy()))
            counter=1
            NS_lists.append(NS_nested_list)
        # Extract the integral factors from the kinetic data for the normalised reaction species with fixed order values. 
        if fixed_order_species!=None:
            counter=0
            for i in fixed_order_species:
                NS_nested_list_fixed_orders=[]
                for exp in experiments:
                    # Obtain the concentration profiles for each of the normalised reaction species i from each experiment exp.
                    NS_nested_list_fixed_orders.append(list(data_scaled[str(exp)][i].to_numpy()))
                NS_lists_fixed_orders.append(NS_nested_list_fixed_orders)
        # Find the number of time points for the kinetic data of each experiment.
        dimensions=[len(l) for l in time_lists]
        # Make the time and normalised species concentration lists from each experiment the same length by padding with 1.
        # Padding is done with 1 rather than 0 to avoid warning messages from trying to raise 0 to negative powers. 
        if all(l==dimensions[0] for l in dimensions)==False:
            max_dim=max(dimensions)
            # Add 0 to the end of each time list which is shorter than the longest list.
            time_lists = [sublist + [1] * (max_dim - len(sublist)) for sublist in time_lists]
            # For each nested list containing the concentration profiles for a given normalised species across all experiments, perform padding. 
            for nl in range(len(NS_lists)):
                NS_lists[nl]=[sublist + [1] * (max_dim - len(sublist)) for sublist in NS_lists[nl]]
            if fixed_order_species!=None:
                # For each nested list containing the concentration profiles for a given normalised species with fixed reaction orders 
                # across all experiments, perform padding. 
                for nl in range(len(NS_lists_fixed_orders)):
                    NS_lists_fixed_orders[nl]=[sublist + [1] * (max_dim - len(sublist)) for sublist in NS_lists_fixed_orders[nl]]
        # Now that the time and normalised species lists are padded, they are converted into 2d and 3d arrays respectively. 
        NS_lists=np.array(NS_lists)
        # Use the convolve function for each 2d array in the 3d array to obtain the average normalised species raised to their respective order values.
        dNS_lists=np.apply_along_axis(lambda x: np.convolve(x, [0.5,0.5], mode='valid'), axis=-1, arr=NS_lists)
        time_lists=np.array(time_lists)
        # Obtain the dt values for each time list.
        dt_lists=np.diff(np.array(time_lists), axis=1)
        # Multiply each dt with the integral fastors from the fixed reaction order reaction species. 
        if fixed_order_species!=None:
            NS_lists_fixed_orders=np.array(NS_lists_fixed_orders)
            # Turn zeros into very small values to avoid errors. 
            dNS_lists_fixed_orders=np.apply_along_axis(lambda x: np.convolve(x, [0.5,0.5], mode='valid'), axis=-1, arr=NS_lists_fixed_orders)
            order_array=np.array(fixed_orders)[:, np.newaxis,np.newaxis]
            # Raise the concentration values of each experiment dataset's normalised species. 
            dNS_lists_fixed_orders_updated=np.power(dNS_lists_fixed_orders, order_array)
            # Multiply together the 2d normalised species arrays and then multiply with the dt values.
            dt_lists=np.apply_along_axis(lambda x: np.prod(x, axis=0), axis=0, arr=dNS_lists_fixed_orders_updated)*dt_lists
        # Access the output concentration data of the first experiment to judge if increasing or decreasing monotonic fitting should be applied. 
        fit_direction='inc'
        if (output_lists[0][0]-output_lists[0][-1])>0:
            fit_direction='dec'
        # Perform first order exploration using "orders". 
        exploration1=explore_orders(dt_lists,output_lists,dNS_lists,dimensions,order_combinations,fit_metric,deg,constraint,fit_direction)
        if exploration1==None:
            return None
        best_orders=exploration1[-1]
        print(f'Best reaction order(s) from the 1st order exploration: {best_orders}')
        # Save the order combination tuples and overlay scores obtained in the first order exploration.
        FQs=exploration1[0]
        order_combinations=exploration1[1]
        # The zoom factor is defined as 140% of the difference between adjacent orders in the first reaction order grid applied if resolution>5. 
        # If the resolution is 4 or 5, the zoom factor is defined as the 110% and 120% of the difference between adjacent orders.  
        mesh_resolution=abs(orders[0]-orders[1])
        multipliers = {x: 1.1 if x == 4 else 1.2 if x == 5 else 1.4 for x in range(4, 10000)}
        zoom_factor=mesh_resolution*multipliers[resolution]
        # Perform a number of additional order explorations given by "interations" to increase the precision of the optimal orders. 
        for i in range(iterations):
            # Define resolution 2 to alternate between odd and even number of values in order value mesh.
            if i//2*2==i:
                resolution2=resolution+1
            else:
                resolution2=resolution
            orders_new=[]
            # Define a new set of order combinations arround the optimal order values from the first order exploration. 
            for k in range(no_orders):
                orders_new.append(np.linspace(best_orders[k]-zoom_factor,best_orders[k]+zoom_factor,resolution2))
            # Perform new order exploration of the new order value combinations. 
            exploration=explore_orders(dt_lists,output_lists,dNS_lists,dimensions,orders_new,fit_metric,deg,constraint,fit_direction)
            if exploration==None:
                return None
            best_orders=exploration[-1]
            print(f'Best reaction order(s) from the {num2words(i+2, to="ordinal_num")} order exploration: {best_orders}')
            # Add datapoints from the latest order exploration to the FQ and order_combinations list. 
            FQs=FQs+exploration[0]
            order_combinations=order_combinations+exploration[1]
            # Define the new zoom factor as before to increase the fineness of the subsequent order value grid. 
            mesh_resolution=abs(exploration[1][-1][-1]-exploration[1][-2][-1])
            zoom_factor=mesh_resolution*multipliers[resolution+1 if resolution2==resolution else resolution]
        # Obtain the minimum point order combination (minium FQ for variance, RMSE and SE, max FQ for R2).
        minFQ=min(FQs)
        if fit_metric=='R2':
            minFQ=max(FQs)
        index_best_order=FQs.index(minFQ)
        best_orders_3=order_combinations[index_best_order]
        # Identify order combinations with overlay scores up to "score_interval" higher than that of the best order:
        if fit_metric=='R2':
            good_order_combos=[order_combinations[i] for i in range(len(order_combinations)) if FQs[i]>minFQ*(1-score_interval)]
        else: 
            good_order_combos=[order_combinations[i] for i in range(len(order_combinations)) if FQs[i]<minFQ*(1+score_interval)]
        # Identify the smallest and largest order values in the good order combinations list for each normalised species.
        largest_values=[max(column) for column in zip(*good_order_combos)]
        smallest_values=[min(column) for column in zip(*good_order_combos)]
        intervals=list(zip(smallest_values, largest_values))
        intervals=np.array(list(zip(smallest_values, largest_values)))
        interval_dict={'Normalised species':normalised_species,f'Lower order limit {score_interval*100}%':intervals[:,0],f'Upper order limit {score_interval*100}%':intervals[:,1]}
        # Create Dataframe for overlay scores against reaction order values.
        result_dict={}
        order_combinations_array=np.array(order_combinations)
        for i in range(len(normalised_species)):
            result_dict[f'order in {normalised_species[i]}']=order_combinations_array[:,i]
        result_dict[f'{fit_metric} overlay score']=FQs
        # Check the final time and print the time taken for the calculation
        end_time=time.time()
        total_time=end_time-start_time
        print(f'The calculation is complete. Time elapsed: {round(total_time,1)} seconds')
        # If the fit degree is set to 0, find the slope of the fit when optimal order values are applied. 
        # The x and y axis are not scaled.
        slope=None
        if deg==1:
            counter=0
            # Update the selection dictonary with the best reaction order values. 
            for species in normalised_species:
                if type(fixed_order_species)==list:
                    if species in fixed_order_species:
                        continue
                selection_copy["normalised_species"][species]=float(best_orders_3[counter])
                counter+=1
            output_species=selection_copy['output_species']
            tables=VTNA_new(data_copy,selection_copy)
            Data_points_y=[]
            Data_points_t=[]
            # Collect all the non-scaled datapoints from the normal VTNA calculation. 
            for i in tables.keys():
                x_values=[float(j) for j in tables[i]['tT'].to_numpy()]
                Data_points_t=Data_points_t+x_values
                ys=[float(j) for j in tables[i][output_species].to_numpy()]
                Data_points_y=Data_points_y+ys
            # Perform the linear fit to determine the slope of the fitted line.  
            if constraint=='origin':
                params = curve_fit(origin_line, Data_points_t, Data_points_y)
                slope=params[0][0]
            else:
                polyfit=np.polyfit(Data_points_t,Data_points_y,deg=deg)
                slope=polyfit[0] 
        return pd.DataFrame(result_dict),best_orders_3,pd.DataFrame(interval_dict),slope

def explore_orders(dt_lists,output_lists,dNS_lists,dimensions,orders,fit_metric,deg,constraint,fit_direction):
    """
    Calculates the normalised time axis of each order combination by calling VTNA_tT() with the \
    condensed kinetic data NS_lists and time_lists. Then performs the specified global fitting \
    procedure to calculate the overlay score as the goodness of fit measure "fit_metric". Returns \
    the overlay scores, order combinations and the order combination giving the best overlay score. 
    Args:
        time_lists (list of lists): Padded list containing time points from each experiment.
        output_lists (list of lists): List containing lists of output concentration profiles from each experiment. 
        NS_lists (list of list of lists): List containing the lists of padded concentration profile lists for \
        each normalised reaction species for each experiment.
        Dimensions (list of int): The number of rows in the DataFrame of each experiment in "data". Used to remove padded values. 
        orders (list of numpy arrays): List of numpy arrays containing the order values to be investigated for each \
        normalised reaction species.
        fit_metric (str): Goodness-of-fit metric chosen to define the overlay score value. Either 'R2', 'RMSE' or 'SE'. \
        otherwise variance is used.
        deg (int): Polynomial degree used in the global fitting prodecure to calculate the overlay score. 
        constraint (str): Defines the constraint of the global fit used by explore_orders() to calculate the overlay score. \
        Can be set to 'monotonic' for monotonic fitting (default), None for ordinary polynomial fitting or 'origin' for \
        linear fit through origin.
        fit_direction (str): The direction of the monotonic fitting procedure, either "inc" or "dec" as determined in \
        VTNA_orders()
    """
    # Define the cancel_calculation as global
    global cancel_calculation
    # Nested function for calculating the normalised time axis
    if constraint=='monotonic':
        # Prepare "polyestimator" for performing monotonic fit with selected degree.
        polyestimator = PolynomRegressor(deg=deg, regularization = None, lam = 0)
        monotone_constraint = Constraints(monotonicity=fit_direction)
    # Define accumulator for goodness of fit measures.
    FQs=[]
    # Get every combination of orders as list of tuples.
    order_pairs=list(itertools.product(*orders))
    # Collapse the output list of lists to use in fitting.
    y=list(itertools.chain(*output_lists))
    for i in order_pairs:
        if cancel_calculation==True:
            cancel_calculation=False
            print("Calculation interrupted by user.")
            return None
        # Obtain the normalised time axis by calling VTNA_tT().
        tT=VTNA_tT(dt_lists,dNS_lists,i).tolist()
        # conncatenate the list of transformed time axis lists while removing padded values using the dimensions list 
        tT=list(itertools.chain(*[tT[i][:dimensions[i]] for i in range(len(dimensions))]))
        # Perform the polynomial fit of "deg" order using the specified fit method and evaluate fitted function at the relevant "tT" values. 
        # Then evaluate the goodness-of-fit by the selected fit metric and add it to the FQs list.
        if constraint=='origin':
            tT_2=[i/max(tT) for i in tT]
            if deg==1:
                params = curve_fit(origin_line, tT_2, y)
                a=params[0]
                y_pred=a*tT_2
                a
            else:
                print('for fit through zero, use deg=1')
                return 
            FQ=score_fit(y,y_pred,fit_metric)
            FQs.append(FQ)
        if constraint=='monotonic': 
            tT_2=np.array([i/max(tT) for i in tT]).reshape((-1,1))    
            polyestimator.fit(tT_2, y, loss = 'l2', constraints={0: monotone_constraint})       
            pred_y = polyestimator.predict(tT_2)
            # Call the score_fit function to evalue the goodness-of-fit using the fit measurement "fit_metric".
            FQ=score_fit(y,pred_y,fit_metric)
            FQs.append(FQ)
        if constraint==None:
            tT_2=[i/max(tT) for i in tT]
            polyfit=np.polyfit(tT_2,y,deg=deg)
            pred_y=np.polyval(polyfit,tT_2)
            FQ=score_fit(y,pred_y,fit_metric)
            FQs.append(FQ)
    # Obtain the best fit metric value (min for SE, RMSE and variance, max for R2)
    minFQ=min(FQs)
    if fit_metric=='R2':
        minFQ=max(FQs)
    # Identify the reaction orderds corresponding to the order minimum point
    index_best_order=FQs.index(minFQ)
    best_orders=order_pairs[index_best_order]
    return FQs, order_pairs, best_orders
    
class CalculationThread(threading.Thread):
# NB: this class function was created and annotated using OpenAI's ChatGTP. 
    def __init__(self, target, args=()):
        # Initialize the thread with a target function and its arguments
        super().__init__(target=target, args=args)
        # Initialize a variable to store the result of the calculation
        self._result = None
    @property
    def result(self):
        # A property method to access the result attribute
        return self._result
    def run(self):
        try:
            # Check if the target function is set
            if self._target:
                # Call the target function with its arguments and store the result
                self._result = self._target(*self._args, **self._kwargs)
        finally:
            # Clean up references to the target function and its arguments
            del self._target, self._args, self._kwargs

def run_calculation_in_thread(window, data, selected_sheet, selected_columns, selected_result_column, calculation_settings):
    global cancel_calculation
    # NB: this function was created and annotated using OpenAI's ChatGTP. 
    # Create a thread for the calculation
    calculation_thread = CalculationThread(target=run_calculation, args=(data, selected_sheet, selected_columns, selected_result_column, calculation_settings))
    calculation_thread.start()
    # Update the GUI to indicate that the calculation is underway
    popup_layout = [[sg.Text("Calculation is underway. Please wait...")],[sg.Button('Stop calculation',key='Cancel')]]
    popup_window = sg.Window("Please Wait", popup_layout, finalize=True)
    # Maintain an infinite loop to process events until the window is shut down.
    cancelled=False
    while True:
        # Read events from the popup window with a timeout of 100 milliseconds
        event, values = popup_window.read(timeout=100)  
        # Check if the popup window is closed or the calculation thread has finished
        if event=='Cancel':
            popup_window.close()
            cancel_calculation=True
            break
        if event == sg.WIN_CLOSED or not calculation_thread.is_alive():
            break
    if cancelled:
        return None
    # Close the popup window
    popup_window.close()
    # Wait for the calculation thread to complete
    calculation_thread.join()
    # Retrieve the return values from the calculation thread
    calculation_result = calculation_thread.result
    return calculation_result

def plot_kinetic_data(data, stochiometry_list, settings):
    '''Generates plots for each experiment in the kinetic dataset "data". 
     If the stochiometry_list argument is specified, plot_data_MB from auto_vtna 
     is used to also show the mass balance for the data.'''
    if settings['mode']=='Together':
        plot_mode='together'
    elif settings['mode']=='Scroll':
        plot_mode='scrollable'
    else:
        plot_mode='separate'
    legend_outside=True if settings['legend_position']=='Outside' else False
    if stochiometry_list != None:
        try:
            # Obtain species and stochiometry lists from dict input.
            species = list(stochiometry_list.keys())
            stochiometry = list(float(i) for i in stochiometry_list.values())
            # plot the data with mass balance. 
            auto_vtna.plot_data_MB(data, species, stochiometry, ylim=settings['ylim'],
                    t_unit=settings['t_unit'],y_unit=settings['y_unit'],
                    fig_size_scaler=settings['size_scaler'],plot_mode=plot_mode,
                    DP_scaler=settings['DP_scaler'],legend_outside=legend_outside,
                    linewidth=settings['linewidth'])
        except ValueError as e:
            sg.popup_error(f"Invalid input: {e}")
    else:
        # Plot the data without mass balance shown. 
        auto_vtna.plot_data(data, ylim=settings['ylim'],
                t_unit=settings['t_unit'],y_unit=settings['y_unit'],
                fig_size_scaler=settings['size_scaler'],plot_mode=plot_mode,
                DP_scaler=settings['DP_scaler'],legend_outside=legend_outside,
                linewidth=settings['linewidth'])

def is_float(value):
    "Checks if a variable can be converted to a float."
    try:
        float(value)
        return True
    except ValueError:
        return False

def is_int(value):
    "Checks if a variable can be converted to an integer."
    try:
        int(value)
        return True
    except ValueError:
        return False

def run_calculation(data, selected_experiments,selected_columns, selected_result_column, calculation_settings):
    '''Runs calculation by preparing the applied calculation settings to be inputted into 
    Auto-VTNA.'''
    # Create a VTNA selection dictonary.
    VTNA_selection_dict=auto_vtna.make_VTNA_selection(data)
    # Add the normalised species selected in the GUI with orders of 1 as placeholders.
    VTNA_selection_dict['normalised_species']={}
    for molecule in selected_columns:
        VTNA_selection_dict['normalised_species'][molecule]=1
    # Define the experiments selected in the GUI that will be included in the calculation.
    VTNA_selection_dict['RM']={}
    for i in selected_experiments:
        VTNA_selection_dict['RM'][i]={}
        VTNA_selection_dict['RM'][i]["omissions"]=None
    # Define the output species specified in the GUI.
    VTNA_selection_dict['output_species']=selected_result_column[0]
    # Define the fixed order species with their selected order values. 
    fixed_order_species=calculation_settings["fixed_order_species"]
    if fixed_order_species!=None:
        for (i,j) in fixed_order_species.items():
            if i in selected_columns:
                VTNA_selection_dict['normalised_species'][i]=j
    # Make a keyword argument for the fixed order species. 
    if fixed_order_species!=None:
        fixed_order_species=list(fixed_order_species.keys())
    # Define the constraint argument based on the selected in teh calculation_settings dictionary. 
    constraint=calculation_settings['constraint']
    if calculation_settings['constraint']=='Via origin':
        constraint='origin'
    if calculation_settings['constraint']=='Monotonic':
        constraint='monotonic'
    if calculation_settings['constraint']=='None':
        constraint=None
    if calculation_settings['initial_mesh_denser']=='False':
        initial_mesh_denser=False
    else:
        initial_mesh_denser=True
    # Set up the automatic VTNA calculation. 
    Calculation_result=Automatic_VTNA(data,VTNA_selection_dict,order_range=calculation_settings["order_range"],
    resolution=calculation_settings["resolution"],deg=calculation_settings["deg"],fit_metric=calculation_settings["fit_metric"],
    iterations=calculation_settings["iterations"], constraint=constraint,score_interval=calculation_settings["score_interval"],
    fixed_order_species=fixed_order_species,initial_mesh_denser=initial_mesh_denser)
    if str(type(Calculation_result.results))!='NoneType':
        VTNA_auto_calculation=Calculation_result
        return VTNA_auto_calculation
    else:
        return None

def check_kinetic_data(kinetic_data):
    '''Code to apply the check_kinetic_data() function from Auto-VTNA 
    and generate a text box for the resulting report.'''
    # NB: this code was written with help from OpenAI's ChatGTP. 
    # Check if the kinetic data is loaded correctly. 
    result_string = auto_vtna.check_kinetic_data(kinetic_data,bold_and_red=False,print_report=False)
    # Create a Tkinter window
    window = tk.Tk()
    window.title("Kinetic Data Check")
    # Create a scrolled text widget to display the result
    text_widget = scrolledtext.ScrolledText(window, wrap=tk.WORD, width=45, height=30)
    text_widget.insert(tk.END, result_string)
    text_widget.configure(state='disabled')  # Make the text read-only
    text_widget.pack(expand=True, fill='both')
    # Run the Tkinter main loop
    window.mainloop()
    return result_string

def add_dash_to_long_words(sentence,length_cutoff):
    '''Function used to add dashes to long words in the initial concentrations
    table. This is required to avoid long column titles.'''
    words = sentence.split()
    for i, word in enumerate(words):
        if len(word) > length_cutoff:
            middle_index = len(word) // 2
            words[i] = word[:middle_index] + "- " + word[middle_index:]
    return ' '.join(words)

def count_significant_figures(number):
    '''Function for finding the number of significant figures, both before and 
    after a '.' to know how whether to apply scientific notation to values in the 
    initial concentration table using apply_best_number_format().'''
    # NB: this function was generated using ChatGTP.
    number_str = str(number).lstrip('0').rstrip('0')
    if '.' in number_str:
        integer_part, decimal_part = number_str.split('.')
        # Count significant figures in the integer part
        integer_figures = len(integer_part)
        # Count significant figures in the decimal part
        decimal_figures = len(decimal_part)
        # Check for leading zeros in the decimal part
        if decimal_part.startswith('0'):
            decimal_figures -= len(decimal_part) - len(decimal_part.lstrip('0'))
        return integer_figures + decimal_figures
    else:
        return len(number_str)
def apply_best_number_format(x, significant_figures):
    '''Checks if the initial concentration value can be rounded or if 
    scientific notation is needed. Then applies the right modification.'''
    number2=round(x,significant_figures)
    significant_figs1=count_significant_figures(x)
    significant_figs2=count_significant_figures(number2)
    if count_significant_figures(number2)<significant_figures and significant_figs1>significant_figs2:
        return f'{x:.{significant_figures}e}'
    else:
        return round(x,significant_figures)

def generate_initial_concentration_table(kinetic_data,fig_size,conc_unit=None,significant_figures=3):
    '''Function to generate an initial concentration popup table based on the 
    information from the initial_concs() function in auto_vtna. The table is made using 
    matplotlib.pyplot'''
    # NB: this code was written with help from OpenAI's ChatGTP.
    # Ensure that all exisiting plots are closed to avoid bugs.
    plt.close('all')
    # Obtain the intial concentrations for each experiment in the kinetic dataset. 
    initial_concs_df= auto_vtna.initial_concs(kinetic_data)
    # Format numeric values with scientific notation and specified significant figures
    initial_concs_df_formatted = initial_concs_df.map(
        lambda x: apply_best_number_format(x, significant_figures) if isinstance(x, (int, float)) else x)
    # Create the pyplot figure.
    fig, ax = plt.subplots(figsize=(12*fig_size, 5*fig_size))
    ax.axis('off')  # Turn off axis
    # Create a table and add it to the axis.
    table = ax.table(cellText=initial_concs_df_formatted.values, colLabels=initial_concs_df_formatted.columns,
     cellLoc='center', loc='center', cellColours=[['#e0e0e0']*len(initial_concs_df_formatted.columns)]*len(initial_concs_df_formatted),
     bbox=[0, 0, 1, 1])
    # Set cell height and width of every cell in the table.
    table.auto_set_column_width(col=list(range(len(initial_concs_df_formatted.columns))))
    for key in table._cells.keys():
        if key[0] == 0:  # First row
            table._cells[key].set_height(0.2)  
        else:
            table._cells[key].set_height(0.1)  
        table._cells[key].set_width(0.015) 
    # Customize appearance for the first row and column to show row and column titles. 
    # Set the cell colour alternate between two shades of grey.
    for i in table._cells:
        # Set the color of every other row following the first to light gray.
        if i[0]!=0 and i[0]%2!=0:
            table._cells[i].set_facecolor('#f8f8f8')
        # Set the color and text of the cells in the first row (column titles)
        if i[0] == 0:  
            table._cells[i].set_facecolor('#6699CC')
            text=table._cells[i].get_text().get_text()
            # Add concentration unit to the column title.
            if conc_unit!=None:
                text=text.replace(' init. conc.','')+f'\nInit.\nconc./{conc_unit}'
            else:
                text=text.replace(' init. conc.','')+f'\nInit.\nconc.'
            # Break down long words in column titles into different lines.
            if len(initial_concs_df.columns)>=6 and len(initial_concs_df.columns)<=10:
                text=add_dash_to_long_words(text,11)
            if len(initial_concs_df.columns)>10:
                text=add_dash_to_long_words(text,7)
            # Use the empty spaces after - to know where to apply line breaks. 
            text=text.replace(' ','\n')
            # Update the text of the cell. 
            table._cells[i].set_text_props(text=text)
        # Set the cells in the first column to blue.
        if i[1]== 0:
            table._cells[i].set_facecolor('#6699CC')
    # Remove the top left corner cell from the table by setting its colour to white.
    table._cells[(0,0)].set_text_props(text='')
    table._cells[(0,0)].set_facecolor('#FFFFFF')
    table._cells[(0,0)].set_linewidth(0)
    # Apply tight layout and show the table. 
    plt.tight_layout()
    plt.show()

def plot_concentration_profiles(kinetic_data, selected_species,settings):
    '''Plots concentration profiles using the plot_data_together function from
    auto_vtna. Applies keyword arguments based on the plot settings dictionary.'''
    legend_outside=True if settings['legend_position']=='Outside' else False
    # Ensure that all exisiting plots are closed to avoid bugs.
    plt.close('all')
    auto_vtna.plot_data_together(kinetic_data, selected_species,ylim=settings['ylim'],
     y_unit=settings['y_unit'], t_unit=settings['t_unit'],legend_outside=legend_outside,
     fig_size_scaler=settings['size_scaler'],DP_scaler=settings['DP_scaler'],linewidth=settings['linewidth'])
    
def main():
    # NB: some parts of the main GUI code was written with help from OpenAI's ChatGTP.
    # Define the cancel_calculation as global
    global cancel_calculation
    img = None # Placeholder for overlay score versus order plot image.
    fixed_orders={} # Placeholder for fixed_orders dictionary
    enable_plot_button = False  # Flag to track whether the calculation has been run
    # Define settings dictionaries for automatic VTNA calculation, the overlay score versus order plot, 
    # data visualisation functions and normal VTNA overlay plot settings. 
    calculation_settings = {"order_range": [-1.5, 2.5], "resolution": 7, "deg": 5, "fit_metric": 'RMSE',
                            "iterations": 7, "constraint": 'Monotonic', "score_interval": 0.15,
                            "fixed_order_species": None,"initial_mesh_denser":'True'}
    order_vs_overlay_plot_settings = {"y_unit":'M',"size_scaler1":1,"popup_scaler":1,"colour_scaler":1,
                            "contour_resolution":None,"datapoints":False,"interval":False,"zoom":'None',
                            "contour_cbar_max":None,"contour_cbar_min":None, "custom_title":None,"show_legend_popup":True,
                            'show_legend_main':True,'decimals_cbar':3,'annotate':True,"score_interval":0.15}
    data_plotting_settings={"ylim":None,"t_unit":None,'y_unit':'M',"size_scaler":1,'DP_scaler':1,"mode":'Scroll',\
                            "SF":3,'legend_position':'Outside','linewidth':1}
    overlay_plot_settings={'y_unit':'M',"size_scaler":1,"tt_scaler":1,"grid":False,"custom_title":None,"check_score":False,\
                           "check_fit":False,"constraint":'Monotonic','deg':5,"fit_metric":'RMSE',
                           'score_settings_visibility':False,'DP_scaler':1,'xtick_rotation':0,'legend':True,
                           'legend_position':'Inside','extra_legend':'True','save':False,'saved_values':{},'line_scaler':1,
                           'tT_notation':'Automatic'}
    # Save the standard settings for resetting the settings window to default values. 
    default_overlay_plot_settings=overlay_plot_settings.copy()
    default_order_vs_overlay_plot_settings=order_vs_overlay_plot_settings.copy()
    default_data_plotting_settings=data_plotting_settings.copy()
    default_calculation_settings=calculation_settings.copy()
    # Define placeholder dict for the Automatic VTNA results. 
    auto_vtna_results={"best_order_dict":{'':'','\t':'','\r':''},"order_vs_overlay_table":None,"interval_table":None,"plot_image":None}
    # Define other miscellaneous placeholders.
    best_order_dict={}
    best_order_dict_rounded={}
    omission_dictionary=None
    omission_range_dictionary={}
    omission_range_dictionary["range_type"]='Percentage'
    omission_mode='Datapoint mode'
    Plot_concentration_profiles_selected_experiments=False
    # Define the layout of the main GUI window.
    layout = [
        [sg.Frame("Load Kinetic Data", layout=[
            [sg.Text("Use Excel file:"),
            sg.InputText(key="file_path", enable_events=True, size=(12, 1)), 
            sg.FileBrowse(),sg.Text("Or:"),sg.Button("Use CSVs",key="Select CSVs")]
        ])],
        [sg.Button("Check data",key="Check Kinetic Data", disabled=True),
         sg.Button("Visualise data",key="Inspect Kinetic Data", disabled=True),
         sg.Button("Crop data", key="Crop Kinetic Data", disabled=True),
         sg.Combo(['Range mode', 'Datapoint mode'], default_value=omission_mode,key="crop_mode", enable_events=True, readonly=True,size=(12,1))],
        [sg.Frame("Select Experiments", layout=[
            [sg.Listbox(values=[], select_mode=sg.LISTBOX_SELECT_MODE_MULTIPLE, size=(22, 5), key="experiments")]
        ]),
        sg.Frame("Select Normalised Species", layout=[
            [sg.Listbox(values=[], select_mode=sg.LISTBOX_SELECT_MODE_MULTIPLE, size=(21, 5), key="columns")]
        ])],
        [sg.Frame("Select Output Species", layout=[
            [sg.Listbox(values=[], select_mode=sg.LISTBOX_SELECT_MODE_SINGLE, size=(20, 1), key="result_column")]
        ]),
        sg.Frame("VTNA Overlay Plot", layout=[
            [sg.Button("Generate plot", key="Generate Overlay Plot", disabled=True),
            sg.Button("Plot settings", key="Overlay Plot Settings", disabled=True)]
        ])],
        [sg.Frame("Automatic VTNA", layout=[
            [sg.Button(" Run ", key='Run Auto VTNA' ,disabled=True),
            sg.Button(" Calculation ⚙", key="calculation_settings", disabled=True),
            sg.Button(" Save \U0001F4BE ", key="Save Results", disabled=not enable_plot_button)],
            [sg.Button("Plot results", disabled=not enable_plot_button),
            sg.Button("Plot ⚙", key="Plot settings", disabled=not enable_plot_button),
            sg.Button("Show +Info", key="extra_info", disabled=not enable_plot_button)]
        ]),
        sg.Table(values=[[key, value] for key, value in auto_vtna_results["best_order_dict"].items()],
                headings=['Compound', 'Order'],
                auto_size_columns=False,
                justification='right',
                num_rows=min(5, len(auto_vtna_results["best_order_dict"])),
                col_widths=[8, 5],
                key='-TABLE-')],
        [sg.Column([
        [sg.pin(sg.Text('placeholder',key="Overlay Score",enable_events=True,visible=False,metadata=False,font=('Any',8)))],
        [sg.pin(sg.Text("Intervals of order values less than",key="Rest of string",enable_events=True,visible=False,metadata=False,font=('Any',8))),
        sg.pin(sg.InputText(key="interval_score", size=(4, 1),default_text=calculation_settings['score_interval']*100,enable_events=True,visible=False,metadata=False,font=('Any',8))),
        sg.pin(sg.Text("% from optimum:",key="Rest of string 2",enable_events=True,visible=False,metadata=False,font=('Any',8)))]], element_justification='left'),
        sg.Column([[sg.pin(sg.Button('Refresh',key='refresh_interval_table',enable_events=True,visible=False,metadata=False))]])],
        [sg.pin(sg.Table(values={},headings=['Normalised species','Lower order limit','Upper order limit'],auto_size_columns=False,col_widths=[14,13,13],key='interval_table',visible=False))]
    ]
    window = sg.Window("Automatic VTNA Calculator", layout)
    # Generate and close a matplotlib pyplot graph to avoid GUI resolution bug.
    plt.plot([1,2],[1,2])
    plt.close()
    # Maintain an infinite loop to process events until the window is shut down
    while True:
        event, values = window.read()
        if event == sg.WIN_CLOSED:
            break
        elif event == "file_path":
            # Update the list of experiments when a file is selected
            file_path = values["file_path"]
            files_calculation=file_path
            if file_path:
                try:
                    # Load the kinetic data as a dictionary of Pandas Dataframe objects. 
                    kinetic_data = pd.read_excel(file_path, sheet_name=None)
                    # Reset the best order dictionary and the saved order values for the overlay plot settings. 
                    best_order_dict_rounded={}
                    overlay_plot_settings['saved_values']={}
                    # Convert column names to strings.
                    for sheet_name, df in kinetic_data.items():
                        df.columns = df.columns.astype(str)
                    # Define placeholders for the data cropping window. 
                    omission_range_dictionary={}
                    omission_range_dictionary["range_type"]='Percentage'
                    # Make a copy of the kinetic data. 
                    data = kinetic_data.copy()
                    # Collect the number of datapoints in each dataframe. 
                    df_nrows=[]
                    for i in kinetic_data.keys():
                        df_nrows.append(kinetic_data[i].shape[0])
                    # Define whether data cropping should operate in datapoint or range mode.
                    # This can be altered by the user. 
                    if max(df_nrows)>20:
                        omission_mode='Range mode'
                    else:
                        omission_mode='Datapoint mode'
                    # Update the "crop_mode" value. 
                    window["crop_mode"].update(value=omission_mode)
                    # Obtain information from the dataset. 
                    time_label=list(data.values())[0].columns[0]
                    experiments = list(data.keys())
                    columns = data[experiments[0]].columns[1:].tolist()  # Exclude the first column
                    # Create a list of reaction species that have concentration profiles that change over time. 
                    reaction_species_not_constant=data[experiments[0]].columns[1:].tolist()
                    for exp in data.keys():
                        for species in columns:
                            # Avoid trying to remove a constant species twice. 
                            if species in reaction_species_not_constant:
                                # Define the list of concentration values for the relevant species in the relevant experiment. 
                                values=list(data[exp][species].to_numpy())
                                # If all values are the same, remove this species from the not constant list
                                if all([species==values[0] for species in values]):
                                    reaction_species_not_constant.remove(species)
                    fixed_order_species = data[experiments[0]].columns[1:].tolist()  # All reaction species
                    # Update the window to add information from the selected kinetic dataset and enable buttons for analysis.
                    window["experiments"].update(values=experiments)
                    window["columns"].update(values=columns)
                    window["result_column"].update(values=reaction_species_not_constant)
                    window["Run Auto VTNA"].update(disabled=False)  # Enable the "Run Automatic VTNA" button
                    window["calculation_settings"].update(disabled=False)  # Enable the "Calculation Settings" button
                    window["Check Kinetic Data"].update(disabled=False)
                    window["Inspect Kinetic Data"].update(disabled=False)
                    window["Generate Overlay Plot"].update(disabled=False)
                    window["Overlay Plot Settings"].update(disabled=False)
                    window["Crop Kinetic Data"].update(disabled=False)
                # Raise any errors:
                except PermissionError as e:
                    sg.popup_error(f"Error opening file (make sure it is not open): {e}")
                    continue
                except ValueError as e:
                    sg.popup_error(f"Error: {e}")
                    continue
                except KeyError as e:
                    sg.popup_error(f"Error: {e}")
                    continue
        # Update the omission_mode variable if the user changes this setting in the main GUI window.
        elif event == "crop_mode":
            omission_mode=values["crop_mode"]
        # Save the results from the automatic VTNA calculation. 
        elif event == "Save Results":
            # Define the layout of the save results pop-up menu. 

            [sg.Column([[sg.Text("Show main plot legend:")],[sg.Text("Show error interval:")],[sg.Text("Overlay score interval: ")],[sg.Text("Main figure size scaler:")]]),
                    sg.Column([[sg.Combo(['True','False'], default_value=str(order_vs_overlay_plot_settings['show_legend_main']),key='show_legend_main', enable_events=True, readonly=True,size=(6,1))],
                               [sg.Combo(['True','False'], default_value=str(order_vs_overlay_plot_settings['interval']),key="interval", enable_events=True, readonly=True,size=(6,1))],
                               [sg.InputText(key="score_interval", size=(6, 1),default_text=order_vs_overlay_plot_settings["score_interval"])],
                               [sg.InputText(key="size_scaler1", size=(6, 1),default_text=float(order_vs_overlay_plot_settings["size_scaler1"]))]])],

            layout_save = [
                [sg.Column([[sg.Text("Filename:")],[sg.Text("File type:")],[sg.Text("Folder location:")]]),
                sg.Column([[sg.InputText(key="filename", size=(5, 1),enable_events=True)],
                [sg.Combo(['csv', 'xlsx'], default_value='xlsx', key="filetype", readonly=True,enable_events=True)],
                [sg.InputText(key="folder", size=(21, 1),enable_events=True,default_text='Default: folder of .exe file')]
                ])],
                [sg.Button('OK'),sg.Button("Select file location")]]
            # Create the save results pop-up menu.
            window_save = sg.Window("Save file settings", layout_save)
            # Apply event loop to process user inputs. 
            while True:
                event_save, values_save = window_save.read()
                if event_save == sg.WINDOW_CLOSED or event_save == 'Cancel':
                    window_save.close()
                    break
                elif event_save == "Select file location":
                    # Create a popup window to allow the user to select file type (csv or xlsx), 
                    # file location and file name. 
                    file_path_save = sg.popup_get_file('Save As', save_as=True, no_window=True,file_types=(("Excel File", "*.xlsx"), ("CSV File", "*.csv")))
                    # Use the inputted file path to define the file name, the extension and the folder location.
                    filename = os.path.splitext(os.path.basename(file_path_save))[0]
                    file_extension = os.path.splitext(file_path_save)[1].replace('.','')
                    folder_location = os.path.dirname(file_path_save)
                    # Update the save window with the relevant info.
                    window_save["filetype"].update(value=file_extension)
                    window_save["filename"].update(value=filename)
                    window_save["folder"].update(value=folder_location)
                elif event_save == 'OK':
                    # use try and exception to enable error messages that don't terminate the GUI.
                    try:
                        if values_save['filename'] is None or values_save['filename'] == '':
                            raise ValueError("Define a filename to continue.")
                        elif values_save['filetype'] is None or window_save['filetype']=='':
                            raise ValueError("Define a filetype to continue.")
                        # Define the file path to use to save the results from the automatic VTNA calculation. 
                        if values_save['folder'] in ['Default: folder of .exe file','']:
                            file_path_save=f"{values_save['filename']}.{values_save['filetype']}"
                        else:
                            # Use the path chosen by the user in the popup_get_file window.
                            file_path_save=f"{values_save['folder']}/{values_save['filename']}.{file_extension}"
                        # Create a report containing information about the automatic VTNA calculation.
                        # First, include the path to the inputted kinetic data file. 
                        files_calculation_nested = [['', files_calculation]]
                        if type(files_calculation) != str:
                            files_calculation_nested = [['', i] for i in files_calculation]
                        files_calculation_nested.insert(0, ["Kinetic data is loaded from:"])
                        # Identify the list of fixed order normalised species. 
                        fixed_order_species_name=VTNA_auto_calculation.fixed_order_species
                        if fixed_order_species_name==None:
                            fixed_order_species_name='None'
                        # Identify the fitting constraint of the calculation. 
                        constraint = VTNA_auto_calculation.constraint if VTNA_auto_calculation.constraint not in [None,''] else 'None'
                        # Create a page containing inputs for the automatic VTNA calculation. 
                        page1 = [
                            ["FITTING METHOD:"],
                            ["  - Fitting constraint:", constraint ],
                            ["  - Degree of fitting polynomial:", VTNA_auto_calculation.deg],
                            ["  - Goodness-of-fit measurement:", VTNA_auto_calculation.fit_metric],
                            ["OTHER CALCULATION SETTINGS:"],
                            ["  - Order value range:", VTNA_auto_calculation.order_range],
                            ["  - Order intervals:", f"{float(VTNA_auto_calculation.score_interval) * 100}% from optimal GOF."],
                            ["  - Algorithm iterations:", VTNA_auto_calculation.iterations],
                            ["  - Density of order value grid:", VTNA_auto_calculation.resolution],
                            ["  - First order grid more dense:",str(VTNA_auto_calculation.initial_mesh_denser)],
                            ["  - Fixed order species:", fixed_order_species_name],
                            ["VTNA SELECTION DICTIONARY:"]
                        ]
                        # Create a list with the contents of the VTNA selection dictionary. 
                        selection_list = list(VTNA_auto_calculation.VTNA_selection.items())
                        # Create a list containing the experiment names and their omissions settings. 
                        s0_2=[['RM']]
                        for key, value in selection_list[0][-1].items():
                            RM_omissions=[key]
                            if list(list(value.items())[0])[1]!=None:
                                RM_omissions=RM_omissions+[list(list(value.items())[0])[0]]+list(list(value.items())[0])[1]
                            else:
                                RM_omissions=RM_omissions+[list(list(value.items())[0])[0]]+['None']
                            s0_2.append(RM_omissions)
                        # Create a list containing all the normalised species defined for the calculation. 
                        s1_2 = [['normalised_species']]
                        s1_2.extend(["  - "+key, value] for key, value in selection_list[2][1].items())
                        selection_list[2] = s1_2
                        selection_list=s0_2+[selection_list[1]]+s1_2
                        # Combine the 3 lists to define the first page of the results file. 
                        page1 = files_calculation_nested + page1 + selection_list
                        # Convert the nested list to a pandas dataframe object. 
                        page1_df = pd.DataFrame(page1)
                    except ValueError as e:
                        sg.popup_error(f"Invalid input: {e}")
                        continue
                    if values_save['filetype'] == 'csv':
                        # Create a report for a CSV file if this has been selected by the user.
                        try:
                            with open(file_path_save, 'w', newline='') as csvfile:
                                csv_writer = csv.writer(csvfile)
                                # Write page1
                                csv_writer.writerows(page1)
                                # Write VTNA_auto_calculation.interval
                                csv_writer.writerow([])
                                csv_writer.writerow([f'{str(float(VTNA_auto_calculation.score_interval)*100)}% order intervals.'])
                                VTNA_auto_calculation.interval.to_csv(csvfile, index=False, header=True)
                                # Write VTNA_auto_calculation.results
                                csv_writer.writerow([])
                                csv_writer.writerow(['Overlay score vs. order matrix'])
                                VTNA_auto_calculation.results.to_csv(csvfile, index=False, header=True)
                        except Exception as e:
                            sg.popup_error(f"Filename permission error (close CSV file with the same name): {e}")
                            continue
                    else:
                        try:
                            # Create a sheet for the page1_df information.
                            with pd.ExcelWriter(file_path_save, engine='xlsxwriter') as writer:
                                page1_df.to_excel(writer, index=False, header=False, sheet_name='Calculation settings')
                            # Create a sheet for the uncertainty intervals of each normalised species.
                            with pd.ExcelWriter(file_path_save, engine='openpyxl', mode='a') as writer:
                                VTNA_auto_calculation.interval.to_excel(writer, index=False, header=True, sheet_name=f'{str(float(VTNA_auto_calculation.score_interval)*100)}% order intervals.')
                            # Create a sheet for the overlay score versus order matrix. 
                            with pd.ExcelWriter(file_path_save, engine='openpyxl', mode='a') as writer:
                                VTNA_auto_calculation.results.to_excel(writer, index=False, header=True, sheet_name='Overlay score vs. order matrix')
                            if img is not None:
                                # Add the image of the recently generated overlay score versus order plot.
                                writer.book.create_sheet('Overlay score vs. order plot')
                                writer.book.save(file_path_save)
                                sheet = writer.sheets['Overlay score vs. order plot']
                                sheet.add_image(img, 'A1')
                                writer.book.save(file_path_save)
                        except Exception as e:
                            sg.popup_error(f"Filename permission error (close excel file with same name): {e}")
                            continue
                    window_save.close()
                    break

        elif event == "Run Auto VTNA":
            # Get selected sheet, columns, and result columns
            selected_sheet = values["experiments"]
            selected_columns = values["columns"]
            selected_result_column = values["result_column"]
            if selected_sheet and selected_columns and selected_result_column:
                # Check that the number of fixed order species is lower than the total number of selected normalised species. 
                count=len(selected_columns)
                if type(calculation_settings["fixed_order_species"])==dict:
                    # Raise an error if any fixed order species are not selected as normalised species. 
                    if not all([i in list(selected_columns) for i in list(calculation_settings["fixed_order_species"].keys())]):
                        sg.popup_error("All selected fixed reaction species must also be selected as normalised species.")
                        continue
                    for i in selected_columns:
                        if i in list(calculation_settings["fixed_order_species"].keys()):
                            count+=-1
                if count>=1:
                        # Run the automatic VTNA calculation in thread to calculate the time passed for the calculation. 
                        calculation = run_calculation_in_thread(window, data, selected_sheet, selected_columns, selected_result_column, calculation_settings)
                        if 'None' not in str(type(calculation.results)):
                            VTNA_auto_calculation = calculation
                            best_order_dict={}
                            # Store the best order information from the calculation. 
                            for i in range(len(VTNA_auto_calculation.best_orders)):
                                best_order_dict[f"{VTNA_auto_calculation.results.columns[i].replace('order in ','')}"]=VTNA_auto_calculation.best_orders[i]
                            orders_last_species=VTNA_auto_calculation.results.iloc[:,-2].to_numpy()
                            resolution_orders=str(format(orders_last_species[-1]-orders_last_species[-2],".12f"))
                            best_order_dict_rounded={}
                            # Make a rounded version of the best order dictionary. 
                            for i in range(len(resolution_orders)):
                                if resolution_orders[i]!='.' and resolution_orders[i]!='0':
                                    break
                            for species,best_order in best_order_dict.items():
                                best_rounded=f'{round(best_order, i-1):.{int(i-1)}f}'
                                best_order_dict_rounded[species]=best_rounded
                            # Update the auto_vtna_result dictionary.
                            auto_vtna_results["best_order_dict"]=best_order_dict
                            auto_vtna_results["order_vs_overlay_table"]=VTNA_auto_calculation.results
                            auto_vtna_results["interval_table"]=VTNA_auto_calculation.interval
                            if VTNA_auto_calculation:
                                # Inform the user that the automatic VTNA calculation is complete. 
                                sg.popup("Calculation complete.",auto_close=True, auto_close_duration=2)
                                window['-TABLE-'].update(values=[[key, value] for key, value in best_order_dict_rounded.items()])
                                enable_plot_button = True  # Enable the "Plot Results" button
                                # Enable buttons for plottin the results of the automatic VTNA calculation. 
                                window["Plot results"].update(disabled=not enable_plot_button)
                                window["Save Results"].update(disabled=not enable_plot_button)
                                window["Plot settings"].update(disabled=not enable_plot_button)
                                window["extra_info"].update(disabled=False)
                                # Identify the best overlay score.
                                if calculation_settings['fit_metric']=='R2':
                                    best_OS=round(max(VTNA_auto_calculation.results[f"{VTNA_auto_calculation.fit_metric} overlay score"]),5)
                                else:
                                    best_OS=round(min(VTNA_auto_calculation.results[f"{VTNA_auto_calculation.fit_metric} overlay score"]),6)  
                                overlay_score_string=f"Overlay score ({calculation_settings['fit_metric']}) at the optimal order values: {str(best_OS)}"
                                if VTNA_auto_calculation.deg==1:
                                    # Determine a sensible number of decimals to include for the slope value. 
                                    rounding=6
                                    slope=VTNA_auto_calculation.best_slope
                                    order_of_magnitude_slope=round(np.log10(slope))
                                    rounding+=(-order_of_magnitude_slope)
                                    # Add the slope to the overlay score string. 
                                    overlay_score_string=overlay_score_string+f'\nSlope of linear fit w/ best orders (w/o data scaling): {round(VTNA_auto_calculation.best_slope,rounding)}'
                                # Update the interval table using the intervals from the automatic VTNA calculation. 
                                df_interval=VTNA_auto_calculation.interval.round(decimals=5)
                                table_values=df_interval.values.tolist()
                                for i,value in enumerate(table_values):
                                    table_values[i]=tuple(value)
                                window['interval_table'].update(values=table_values)
                                window["Overlay Score"].update(overlay_score_string)
                                window["interval_score"].update(calculation_settings['score_interval']*100)
                                window['interval_table'].update(num_rows=len(table_values))   
                else:
                    sg.popup_error("The number of fixed order species can't be the same as or higher than the number of normalised reaction species.")
            else:
                sg.popup_error("Please select experiments, reaction species, and output reaction species to run automatic VTNA.")
        elif event == "refresh_interval_table":
            try:
                interval_score=values['interval_score']
                if not is_float(interval_score):
                    raise ValueError("Invalid score interval. Must be numerical")
                if not float(interval_score)>0:
                    raise ValueError("Invalid score interval. Must be above 0.")
            except ValueError as e:
                sg.popup_error(f"Invalid input: {e}")
                continue
            intervals=VTNA_auto_calculation.identify_error_intervals(score_interval=float(float(interval_score)/100)).round(decimals=5)
            table_values=intervals.values.tolist()
            for i,value in enumerate(table_values):
                table_values[i]=tuple(value)
            window['interval_table'].update(values=table_values)
        elif event == 'extra_info':
            # Update the extra info GUI window extension status.
            element_table = window['interval_table']
            # Switch the visibility status of the element table. 
            visible = element_table.metadata = not element_table.metadata
            window['interval_table'].update(visible=visible)
            window["Overlay Score"].update(visible=visible)
            window["Rest of string"].update(visible=visible)
            window["Rest of string 2"].update(visible=visible)
            window["interval_score"].update(visible=visible)
            window['refresh_interval_table'].update(visible=visible)
            # Update the extra info button according to the visibility status. 
            if visible:
                window['extra_info'].update('Hide + Info')
            else:
                window['extra_info'].update('Show +Info')

        elif event == "Plot settings":
            # Define the layout of the plot settings menu, consisting of 3 different frames for general settings, 
            # contour plot settings and pop-up overlay plot settings (for the VTNA overlay plots generated by clicking 
            # on the main plot.)
            layout_settings_OvsO_plot = [
                [sg.Text("Overlay Score vs. Order Plot Settings", font=("helvetica", 11, "bold"))],
                [sg.Frame("General settings",layout=[
                    [sg.Column([[sg.Text("Show main plot legend:")],[sg.Text("Show error interval:")],[sg.Text("Overlay score interval: ")],[sg.Text("Main figure size scaler:")]]),
                    sg.Column([[sg.Combo(['True','False'], default_value=str(order_vs_overlay_plot_settings['show_legend_main']),key='show_legend_main', enable_events=True, readonly=True,size=(6,1))],
                               [sg.Combo(['True','False'], default_value=str(order_vs_overlay_plot_settings['interval']),key="interval", enable_events=True, readonly=True,size=(6,1))],
                               [sg.InputText(key="score_interval", size=(6, 1),default_text=order_vs_overlay_plot_settings["score_interval"])],
                               [sg.InputText(key="size_scaler1", size=(6, 1),default_text=float(order_vs_overlay_plot_settings["size_scaler1"]))]])],
                    [sg.Text("Overlay score range:"),sg.InputText(key="contour_cbar_min", size=(6, 1),default_text='Default' if order_vs_overlay_plot_settings["contour_cbar_min"] is None else str(order_vs_overlay_plot_settings["contour_cbar_min"])),
                    sg.Text("to"),sg.InputText(key="contour_cbar_max", size=(6, 1),default_text='Default' if order_vs_overlay_plot_settings["contour_cbar_max"] is None else str(order_vs_overlay_plot_settings["contour_cbar_max"]))],
                    [sg.Text("Custom figure title:"),
                    sg.InputText(key="custom_title", size=(18, 1),default_text=str(order_vs_overlay_plot_settings["custom_title"]) if order_vs_overlay_plot_settings["custom_title"]!=None else 'Default')],
                    ])],
                [sg.Frame("Settings for contour plot",layout=[
                    [sg.Column([[sg.Text("Show datapoints:")],[sg.Text("Show optimum annotation:")],[sg.Text("Color bar max scaler:")],[sg.Text("Colour bar resolution:")],[sg.Text("Zoom around optimum:")],[sg.Text("Number of cbar decimals:")]]),
                    sg.Column([[sg.Combo(['True','False'], default_value=str(order_vs_overlay_plot_settings['datapoints']),key="datapoints", enable_events=True, readonly=True,size=(6,1))],
                               [sg.Combo(['True','False'], default_value=str(order_vs_overlay_plot_settings['annotate']),key="annotate", enable_events=True, readonly=True,size=(6,1))],
                               [sg.InputText(key="colour_scaler", size=(6, 1),default_text=str(order_vs_overlay_plot_settings["colour_scaler"]), enable_events=True)],
                               [sg.InputText(key="contour_resolution", size=(6, 1),default_text='Default' if order_vs_overlay_plot_settings["contour_resolution"] is None else str(order_vs_overlay_plot_settings["contour_resolution"]), enable_events=True)],
                               [sg.InputText(key="zoom", size=(6, 1),default_text=str(order_vs_overlay_plot_settings["zoom"]), enable_events=True)],
                               [sg.InputText(key="decimals_cbar", size=(6, 1),default_text=str(order_vs_overlay_plot_settings["decimals_cbar"]), enable_events=True)]])]])],
                [sg.Frame("Settings for pop-up overlay plots",layout=[
                    [sg.Column([[sg.Text("Show legend:")],[sg.Text("Pop-up figure size scaler:")],[sg.Text("Concentration unit:")]]),
                    sg.Column([[sg.Combo(['True','False'], default_value=str(order_vs_overlay_plot_settings['show_legend_popup']),key="show_legend_popup", enable_events=True, readonly=True,size=(6,1))],
                               [sg.InputText(key="popup_scaler", size=(6, 1),default_text=str(order_vs_overlay_plot_settings["popup_scaler"]),enable_events=True)],
                               [sg.InputText(key="y_unit", size=(6, 1),default_text=str(order_vs_overlay_plot_settings["y_unit"]), enable_events=True)]])]])],
                [sg.Button("OK"), sg.Button("Cancel"),sg.Button("Reset settings")]]
            # Create the plot settings window. 
            window_OvsO_plot_settings = sg.Window("Order Versus Overlay Plot Settings", layout_settings_OvsO_plot)
            # Maintain an infinite loop to process events until the window is shut down.
            while True:
                event_OvsO_plot_settings, values_OvsO_plot_settings = window_OvsO_plot_settings.read()
                # Break the while loop to close the window if cancel is pressed or the window is closed. 
                if event_OvsO_plot_settings == sg.WIN_CLOSED or event_OvsO_plot_settings == "Cancel":
                    break
                # If the reset button is pressed, reset all the settings.
                if event_OvsO_plot_settings=="Reset settings":
                    for i,j in default_order_vs_overlay_plot_settings.items():
                        if j==None:
                            j='Default'
                        elif 'scaler' in str(i):
                            j=float(j)
                        else:
                            j=str(j)
                        window_OvsO_plot_settings[i].update(value=j)
                # if the user clicks OK, save the settings selected. 
                elif event_OvsO_plot_settings == "OK":
                    try:
                        order_vs_overlay_plot_settings["y_unit"] = str(values_OvsO_plot_settings["y_unit"])
                        main_size_scaler = float(values_OvsO_plot_settings["size_scaler1"])
                        popup_scaler = float(values_OvsO_plot_settings["popup_scaler"])
                        colour_scaler = float(values_OvsO_plot_settings["colour_scaler"])
                        if values_OvsO_plot_settings["contour_resolution"] not in [None,'Default']:
                            order_vs_overlay_plot_settings["contour_resolution"]=float(values_OvsO_plot_settings["contour_resolution"])
                        else:
                            order_vs_overlay_plot_settings["contour_resolution"]=None 
                        datapoints = values_OvsO_plot_settings["datapoints"]
                        annotate = values_OvsO_plot_settings["annotate"]
                        interval = values_OvsO_plot_settings["interval"]
                        plot_score_interval=values_OvsO_plot_settings["score_interval"]
                        show_legend_popup = values_OvsO_plot_settings["show_legend_popup"]
                        show_legend_main = values_OvsO_plot_settings["show_legend_main"]
                        zoom_range = str(values_OvsO_plot_settings["zoom"])
                        decimals_cbar = values_OvsO_plot_settings["decimals_cbar"]
                        order_vs_overlay_plot_settings["decimals_cbar"]=decimals_cbar
                        if values_OvsO_plot_settings["contour_cbar_max"] not in [None,'Default','']:
                            order_vs_overlay_plot_settings["contour_cbar_max"] = float(values_OvsO_plot_settings["contour_cbar_max"])
                        else:
                            order_vs_overlay_plot_settings["contour_cbar_max"] = None
                        if values_OvsO_plot_settings["contour_cbar_min"] not in [None,'Default','']:
                            order_vs_overlay_plot_settings["contour_cbar_min"] = float(values_OvsO_plot_settings["contour_cbar_min"])
                        else:
                            order_vs_overlay_plot_settings["contour_cbar_min"]=None
                        # Define the datapoint and interval settings by bolean values.
                        if datapoints=='True':
                            order_vs_overlay_plot_settings["datapoints"] = True
                        else:
                            order_vs_overlay_plot_settings["datapoints"] = False
                        if annotate=='True':
                            order_vs_overlay_plot_settings["annotate"] = True
                        else:
                            order_vs_overlay_plot_settings["annotate"] = False
                        if interval=='True':
                            order_vs_overlay_plot_settings["interval"] = True
                        else:
                            order_vs_overlay_plot_settings["interval"] = False
                        if show_legend_popup=='True':
                            order_vs_overlay_plot_settings["show_legend_popup"] = True
                        else: 
                            order_vs_overlay_plot_settings["show_legend_popup"] = False
                        if show_legend_main=='True':
                            order_vs_overlay_plot_settings["show_legend_main"] = True
                        else: 
                            order_vs_overlay_plot_settings["show_legend_main"] = False
                        # Check for erroneous inputs. 
                        if not float(plot_score_interval)>0:
                            raise ValueError("Invalid score interval for plot. Must be above 0.")
                        if not (0.1 < main_size_scaler <= 3):
                            raise ValueError("Invalid figure size scaler. Must be between 0.1 and 3.")
                        if not (0.1 < popup_scaler <= 3):
                            raise ValueError("Invalid popup figure size scaler. Must be between 0.1 and 3.")
                        if not (0 < colour_scaler <= 2):
                            raise ValueError("Invalid contour plot colour bar scaler. Must be between 0 and 2")
                        if float(plot_score_interval)<=0:
                            raise ValueError("Invalid overlay score interval. Must be above 0.")
                        if zoom_range!='None':
                            if not (0 < float(zoom_range) <= float(calculation_settings["order_range"][-1]-calculation_settings["order_range"][0])+2):
                                raise ValueError("The zoom range defines the breadth of order values shown in the\
 order versus overlay plot, and can't be below 0 and should't be much higher than the order range of the calculation")
                        # Save the relevant settings if no error was produced. 
                        order_vs_overlay_plot_settings["size_scaler1"] = main_size_scaler
                        order_vs_overlay_plot_settings["popup_scaler"] = popup_scaler
                        order_vs_overlay_plot_settings["colour_scaler"] = colour_scaler
                        order_vs_overlay_plot_settings["score_interval"]=float(plot_score_interval)
                        if zoom_range=='None' or len(zoom_range)==0:
                            order_vs_overlay_plot_settings["zoom"]='None'
                        else:
                            order_vs_overlay_plot_settings["zoom"]= float(zoom_range)
                    except ValueError as e:
                        sg.popup_error(f"Invalid input: {e}")
                        continue
                    window_OvsO_plot_settings.close()
            window_OvsO_plot_settings.close()

        elif event == "Overlay Plot Settings":
            # Create the layout for the plot settings menu. 
            # The pinned settings are for when the user selects to show the fitted curve and or overlay score value. 
            # This opens up extra settings like which polynomial degree and which constraint to apply. 
            layout_settings = [
                [sg.Text("Overlay Plot Settings", font=("helvetica", 12, "bold"))],
                [sg.Frame("Axis Settings", layout=[
                [sg.Text("Conc. unit:"),
                sg.InputText(key="y_unit", size=(4, 1), default_text=str(overlay_plot_settings["y_unit"])),
                sg.Text("x-tick label rotation:"),
                sg.InputText(key="xtick_rotation", size=(4, 1), default_text=float(overlay_plot_settings["xtick_rotation"]))
                ],
                [sg.Text("x-axis number type:"),
                sg.Combo(['Automatic','Normal', 'Scientific'],size=(9, 1), default_value=str(overlay_plot_settings['tT_notation']), key="tT_notation", enable_events=True, readonly=True)]
                ])],
                [sg.Frame("Size Scalers", layout=[
                [sg.Text("Linewidth scaler:"),
                sg.InputText(key="line_scaler", size=(4, 1), default_text=float(overlay_plot_settings["line_scaler"])),
                sg.Text("Figure scaler:"),
                sg.InputText(key="size_scaler", size=(4, 1), default_text=float(overlay_plot_settings["size_scaler"]))],
                [sg.Text("Datapoint scaler:"),
                sg.InputText(key="DP_scaler", size=(4, 1), default_text=float(overlay_plot_settings["DP_scaler"])),
                sg.Text("x-label scaler:"),
                sg.InputText(key="tt_scaler", size=(4, 1), default_text=float(overlay_plot_settings["tt_scaler"]))]
                ])],

                [sg.Frame("Other Settings", layout=[
                [sg.Text("Show legend: "),
                sg.Combo(['True', 'False'],size=(5, 1), default_value=str(overlay_plot_settings['legend']), key="legend", enable_events=True, readonly=True),
                sg.Text("Gridlines:"),
                sg.Combo(['True', 'False'],size=(5, 1), default_value=str(overlay_plot_settings['grid']), key="grid", enable_events=True, readonly=True)],
                [sg.Text("Legend position:"),
                sg.Combo(['Inside','Outside'],size=(5, 1), default_value=str(overlay_plot_settings['legend_position']), key="legend_position", enable_events=True, readonly=True)],
                [sg.Text("Custom title:"),
                sg.InputText(key="custom_title", size=(24, 1),default_text=str(overlay_plot_settings["custom_title"]) if overlay_plot_settings["custom_title"]!=None else 'Default')]
                ])],
                [sg.Checkbox("Calculate overlay score", key='check_score',default=overlay_plot_settings["check_score"],enable_events=True),
                    sg.Checkbox("Show fit function", key='check_fit',default=overlay_plot_settings["check_fit"],enable_events=True)],
                [sg.pin(sg.Text("Constraint:",enable_events=True,visible=overlay_plot_settings['score_settings_visibility'],key="constraint_text")),sg.pin(sg.Combo(['Monotonic', 'None', 'Via origin'],
                        default_value=overlay_plot_settings["constraint"],size=(9, 1),
                        key="constraint", enable_events=True, readonly=True,visible=overlay_plot_settings['score_settings_visibility'])),
                sg.pin(sg.Text("Degree:",enable_events=True,visible=overlay_plot_settings['score_settings_visibility'],key='deg_text')),sg.pin(sg.Combo([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], default_value=overlay_plot_settings["deg"],
                        key="deg",enable_events=True, readonly=True,visible=overlay_plot_settings['score_settings_visibility']))],
                [sg.pin(sg.Text("Goodness-of-fit measurement:",enable_events=True,visible=overlay_plot_settings['score_settings_visibility'],key='GOF_text')),sg.pin(sg.Combo(['SE', 'RMSE', 'R2', 'Variance'],
                        default_value=overlay_plot_settings["fit_metric"],size=(8, 1),
                        key="fit_metric", enable_events=True, readonly=True,visible=overlay_plot_settings['score_settings_visibility']))],
                [sg.pin(sg.Text("Show fit function legend:",enable_events=True,visible=overlay_plot_settings['score_settings_visibility'],key='extra_legend_text')),sg.pin(sg.Combo(['True','False'], default_value=overlay_plot_settings["extra_legend"],
                        key="extra_legend",enable_events=True, readonly=True,visible=overlay_plot_settings['score_settings_visibility']))],
                [sg.Button("OK"), sg.Button("Cancel"), sg.Button("Reset settings")]
            ]
            # Create the overlay plot settings window. 
            window_overlay_settings = sg.Window("Overlay Plot Settings", layout_settings,resizable=True)
            # Maintain an infinite loop to process events until the window is shut down.
            while True:
                event_overlay_settings, values_overlay_settings = window_overlay_settings.read()
                if event_overlay_settings == sg.WIN_CLOSED or event_overlay_settings == "Cancel":
                    break
                # Update to the default settings if the user presses the reset settings button. 
                elif event_overlay_settings=='Reset settings':
                    # Loop through the elements of the default overlay plot settings and apply these
                    # to the overlay settings window.
                    for i,j in default_overlay_plot_settings.items():
                        if 'scaler' in i:
                            j=float(j)
                        elif i=='custom_title':
                            j='Default'
                        else:
                            j=str(j)
                        if i not in ['score_settings_visibility','save','saved_values']:
                            window_overlay_settings[i].update(value=j)
                    window_overlay_settings['check_fit'].update(value=False)
                    window_overlay_settings['check_score'].update(value=False)
                    # Set all the extra settings to be invisible
                    overlay_plot_settings['score_settings_visibility']=False
                    window_overlay_settings["constraint_text"].update(visible=False)
                    window_overlay_settings["constraint"].update(visible=False)
                    window_overlay_settings["deg_text"].update(visible=False)
                    window_overlay_settings["deg"].update(visible=False)
                    window_overlay_settings["GOF_text"].update(visible=False)
                    window_overlay_settings["fit_metric"].update(visible=False)
                    window_overlay_settings["extra_legend"].update(visible=False)
                    window_overlay_settings["extra_legend_text"].update(visible=False)
                if event_overlay_settings in ['check_score','check_fit']:
                    # If one of the check boxes for overlay score or show the fit curve have been modified, check whether the
                    # boxes are ticked (values True). If one or both of these are True, define the visibility variable as True and 
                    # use this to make the extra settings visible. 
                    if values_overlay_settings['check_score']==True or values_overlay_settings['check_fit']==True:
                        visibility=True
                    else:
                        visibility=False
                    overlay_plot_settings['score_settings_visibility']=visibility
                    window_overlay_settings["constraint_text"].update(visible=visibility)
                    window_overlay_settings["constraint"].update(visible=visibility)
                    window_overlay_settings["deg_text"].update(visible=visibility)
                    window_overlay_settings["deg"].update(visible=visibility)
                    window_overlay_settings["GOF_text"].update(visible=visibility)
                    window_overlay_settings["fit_metric"].update(visible=visibility)
                    window_overlay_settings["extra_legend"].update(visible=visibility)
                    window_overlay_settings["extra_legend_text"].update(visible=visibility)
                elif event_overlay_settings == "OK":
                    # Try to define the settings variables based on the overlay settings window values. 
                    try:
                        y_unit = values_overlay_settings["y_unit"]
                        tt_scaler = float(values_overlay_settings["tt_scaler"])
                        DP_scaler = float(values_overlay_settings["DP_scaler"])
                        line_scaler = float(values_overlay_settings["line_scaler"])
                        xtick_rotation = float(values_overlay_settings["xtick_rotation"])
                        size_scaler = float(values_overlay_settings["size_scaler"])
                        grid = values_overlay_settings["grid"]
                        legend_bol = values_overlay_settings["legend"]
                        overlay_plot_settings["legend_position"] = values_overlay_settings["legend_position"]
                        overlay_plot_settings["custom_title"]=values_overlay_settings["custom_title"]
                        overlay_plot_settings["size_scaler"]=size_scaler
                        overlay_plot_settings["y_unit"]=y_unit
                        overlay_plot_settings["tt_scaler"]=tt_scaler
                        overlay_plot_settings["DP_scaler"]=DP_scaler
                        overlay_plot_settings["line_scaler"]=line_scaler
                        deg=values_overlay_settings["deg"]
                        constraint=values_overlay_settings["constraint"]
                        overlay_plot_settings["fit_metric"]=values_overlay_settings["fit_metric"]
                        overlay_plot_settings["check_score"]=values_overlay_settings["check_score"]
                        overlay_plot_settings["check_fit"]=values_overlay_settings["check_fit"]
                        overlay_plot_settings["tT_notation"]=values_overlay_settings["tT_notation"]
                        # Define string boleans as proper boleans. 
                        if grid=='True':
                            overlay_plot_settings["grid"]=True
                        else:
                            overlay_plot_settings["grid"]=False
                        if legend_bol=='True':
                            overlay_plot_settings["legend"]=True
                        else:
                            overlay_plot_settings["legend"]=False
                        # Check if the selected values are valid. Raise errors if not. 
                        if deg not in [1, 2, 3, 4, 5, 6, 7,8,9,10,11,12,13]:
                            raise ValueError("Invalid degree. Must be an integer between 1 and 13.")
                        if constraint=='Via origin':
                            if int(deg)!=1:
                                raise ValueError("Set fitting degree to 1 for linear line via origin.")
                        if not (0.1 < size_scaler <= 3):
                            raise ValueError("Invalid figure size scaler. Must be between 0.1 and 3.")
                        if not (0.1 < tt_scaler <= 3):
                            raise ValueError("Invalid figure transformed time axis label scaler. Must be between 0.1 and 3.")
                        if not (0.00 <= DP_scaler <= 5):
                            raise ValueError("Invalid datapoint size scaler. Must be between 0 and 5.")
                        if not (0.00 <= xtick_rotation <= 180):
                            raise ValueError("Invalid x-axis tick label rotation value. Must be between 0 and 180.")
                        # Update the overlay plot settings now that they have been checked for errors. 
                        overlay_plot_settings["deg"]=deg
                        overlay_plot_settings["constraint"]=constraint
                        overlay_plot_settings["extra_legend"]=values_overlay_settings["extra_legend"]
                        overlay_plot_settings["xtick_rotation"]=xtick_rotation
                    except ValueError as e:
                            sg.popup_error(f"Invalid input: {e}")
                            continue
                    window_overlay_settings.close()
            window_overlay_settings.close()

        elif event == "Plot results":
            # Plot the results from the automatic VTNA calculation using settings from the over_vs_overlay_plot_settings dictionary. 
            try:
                # Check the dimensionality of the overlay score versus order matrix to see if a graph (one dimension) or contour plot 
                # (two dimensions) can be produced.
                no_dimensions=-len(fixed_orders.values())+len(selected_columns)
                if no_dimensions>2:
                    raise ValueError("Can't plot results for more than two normalised reaction species to investigate.")
                # Ensure that all exisiting plots are closed to avoid bugs.
                plt.close('all')
                # Define the zoom range and custom title inputs. 
                zoom_range=order_vs_overlay_plot_settings["zoom"]
                if zoom_range=='None':
                    zoom_range=False
                custom_title = order_vs_overlay_plot_settings["custom_title"]
                if custom_title!=None and custom_title!='Default':
                    if len(custom_title)==0 or custom_title=='None':
                        custom_title=None
                # Create the figure and save it as a variable img.
                img=VTNA_auto_calculation.plot_orders_vs_overlay(y_unit=order_vs_overlay_plot_settings["y_unit"],size_scaler=order_vs_overlay_plot_settings["size_scaler1"],
                size_scaler2=order_vs_overlay_plot_settings["popup_scaler"],color_scaler=order_vs_overlay_plot_settings["colour_scaler"],
                fixed_cbar_resolution=order_vs_overlay_plot_settings["contour_resolution"],points=order_vs_overlay_plot_settings["datapoints"],
                interval=order_vs_overlay_plot_settings["interval"],zoom=zoom_range,overlay_score_range_max=order_vs_overlay_plot_settings["contour_cbar_max"],
                overlay_score_range_min=order_vs_overlay_plot_settings["contour_cbar_min"],title=custom_title,show_legend_popup=order_vs_overlay_plot_settings["show_legend_popup"],
                show_legend_main=order_vs_overlay_plot_settings["show_legend_main"],decimals_in_colorbar=order_vs_overlay_plot_settings["decimals_cbar"],
                annotate_optimum=order_vs_overlay_plot_settings["annotate"],specified_score_interval=order_vs_overlay_plot_settings['score_interval'])
            except ValueError as e:
                sg.popup_error(f"Invalid input: {e}")
                continue
        
        elif event == "Generate Overlay Plot":
            # Ensure that all exisiting plots are closed to avoid bugs.
            plt.close('all')
            # Define the selected experiments, selected columns and selected output species
            # for generating the correct VTNA overlay plot. 
            selected_sheet = values["experiments"]
            selected_columns = values["columns"]
            selected_result_column = values["result_column"]
            # Define keyword arguments from the overlay_plot_settings dictionary. 
            tt_scaler=overlay_plot_settings['tt_scaler']
            DP_scaler=float(overlay_plot_settings["DP_scaler"])
            line_scaler=float(overlay_plot_settings["line_scaler"])
            xtick_rotation=float(overlay_plot_settings["xtick_rotation"])
            y_unit=overlay_plot_settings['y_unit']
            size_scaler=overlay_plot_settings['size_scaler']
            grid=overlay_plot_settings['grid']
            tT_notation=overlay_plot_settings['tT_notation']
            title=None if overlay_plot_settings['custom_title']=='Default' else overlay_plot_settings['custom_title']
            constraint=overlay_plot_settings['constraint']
            extra_legend=True if overlay_plot_settings['extra_legend']=='True' else False
            fit_metric=overlay_plot_settings['fit_metric']
            deg=overlay_plot_settings['deg']
            show_overlay_score=overlay_plot_settings['check_score']
            show_fit_function=overlay_plot_settings['check_fit']
            legend_bol=overlay_plot_settings['legend']
            if overlay_plot_settings['legend_position']=='Inside':
                legend_outside=False
            else:
                legend_outside=True
            if overlay_plot_settings['custom_title'] in ['','None']:
                title=None
            if selected_sheet and selected_columns and selected_result_column:
                # Open a pop-up window to input order values for selected fixed order species
                if selected_columns:
                    layout_order_values = [
                        [sg.Text(f"Enter order values for selected species:")],
                    ]
                    # Add a text and input text element for each normalised species to allow the user to input order values. 
                    species_column=[]
                    orders_column=[]
                    for species in selected_columns:
                        species_column.append([sg.Text(f"{species}:")])
                        if overlay_plot_settings['save'] and species in list(overlay_plot_settings['saved_values'].keys()):
                            orders_column.append([sg.InputText(key=f"order_value_{species}", size=(5, 1),enable_events=True,default_text=overlay_plot_settings['saved_values'][species])])
                        else:
                            orders_column.append([sg.InputText(key=f"order_value_{species}", size=(5, 1),enable_events=True)])
                    layout_order_values.append([sg.Column(species_column,element_justification='left'),sg.Column(orders_column,element_justification='left')])
                    # Define the OK and cancel buttons for the pop-up menu.
                    buttons=[sg.Button("OK"), sg.Button("Cancel")]
                    # If optimal order values have been calculated by automatic VTNA, and some of these reaction species have 
                    # been selected for the overlay plot, create a button for auto-filling these calculated order values. 
                    if len(best_order_dict)!=0:
                        if [i in selected_columns for i in list(best_order_dict.keys())]:
                            buttons.append(sg.Button("Auto-fill calculated values.",key='auto'))
                    layout_order_values.append(buttons)
                    # Add a save order values checkbox that will allow the user to keep the inputted order values even after 
                    # closing the pop-up menu.  
                    layout_order_values.append([sg.Text("Save order values:"),sg.Checkbox("", key="save", enable_events=True,
                        default=overlay_plot_settings['save']),sg.Button("Reset orders")])
                    window_order_values = sg.Window("Enter Order Values", layout_order_values,resizable=True,finalize=True)
                    # Maintain an infinite loop to process events until the window is shut down.
                    while True:
                        input_valid = False
                        event_order_values, values_order_values = window_order_values.read()
                        if event_order_values == sg.WIN_CLOSED or event_order_values == "Cancel":
                            # User closed the window or pressed "Cancel"
                            plt.close('all')
                            window_order_values.close()
                            break
                        if 'order_value_' in event_order_values:
                            # Update the saved values dictionary in the overlay_plot_settings dictionary if a value is altered.
                            if overlay_plot_settings['save']:
                                order_values={}
                                for species in selected_columns:
                                    order_values[species] = str(values_order_values[f"order_value_{species}"])
                                overlay_plot_settings['saved_values']=order_values
                        # Update the "save" setting if the user clicks the tick box.
                        if event_order_values=='save':
                            overlay_plot_settings['save']=not overlay_plot_settings['save']
                            # If the save setting is set to True, save the inputted order values to the overlay_plot_settings dictionary. 
                            if overlay_plot_settings['save']:
                                order_values={}
                                for species in selected_columns:
                                    order_values[species] = str(values_order_values[f"order_value_{species}"])
                                overlay_plot_settings['saved_values']=order_values
                        # If the user clicks the auto fill button, update the relevant order value input buttons. 
                        if event_order_values=='auto':
                            for species,order in best_order_dict_rounded.items():
                                window_order_values[f"order_value_{species}"].update(value=order)
                        # If the user clicks the reset orders button, update all order values to be empty. 
                        if event_order_values=="Reset orders":
                            for species in selected_columns:
                                window_order_values[f"order_value_{species}"].update(value="")
                        if event_order_values=="OK":
                            order_values = {}
                            # Check that the inputted order value is valid for each of the normalised reaction species. 
                            for species in selected_columns:
                                input_valid = True
                                try:
                                    order_value = float(values_order_values[f"order_value_{species}"])
                                    # Check if the value is within the specified range
                                    if not (-4 <= order_value <= 4):
                                        sg.popup_error(f"Invalid input for order value of {species}. Must be between -4 and 4.")
                                        input_valid = False
                                        break
                                except ValueError:
                                    sg.popup_error(f"Invalid input for order value of {species}. Must be a numerical value.")
                                    input_valid = False
                                    break
                                order_values[species] = str(values_order_values[f"order_value_{species}"])
                            # If input is valid, break out of the loop
                            if input_valid:
                                # Save the final order values to the overlay_plot_settings dictionary. 
                                if overlay_plot_settings['save']:
                                    overlay_plot_settings['saved_values']=order_values
                                # Generate the overlay plot using the overlay_plot function. 
                                overlay_plot(data,selected_sheet,order_values,selected_result_column[0],tt_scaler=tt_scaler,\
                                         grid=grid,y_unit=y_unit,size_scaler=size_scaler,title=title, fit_metric=fit_metric,
                                         deg=deg,constraint=constraint,show_overlay_score=show_overlay_score,\
                                         show_fit_function=show_fit_function,DP_scaler=DP_scaler,xtick_rotation=xtick_rotation,\
                                         show_legend=legend_bol,legend_outside=legend_outside,extra_legend=extra_legend,\
                                         line_scaler=line_scaler,tT_notation=tT_notation)
            else:
                sg.popup_error("Please select experiments, reaction species, and output reaction species before generating overlay plot.")

        elif event == "calculation_settings":
            # Create the layout for fixed order species based on the pre-existing fixed order species 
            # dictionary in the calculation settings dictionary. 
            fixed_orders={}
            layout_fixed_order=[]
            label_width = max(len(str(compound)) for compound in fixed_order_species)
            fixed_orders_start=calculation_settings["fixed_order_species"]
            if type(fixed_orders_start)!=dict:
                fixed_orders_start={}
            # Loop through each compound column in the dataset and make an entry in the fixed order species frame.
            for compound in fixed_order_species:
                if compound in fixed_orders_start.keys():
                        # If the compound is already present in the fixed order species dictionary, 
                        # check the checkbox and fill in the order values. 
                        layout_fixed_order.append([
                            sg.Text(compound, size=(label_width, 1)),
                            sg.Checkbox("", key=f"Check_{compound}", enable_events=True,default=True),
                            sg.InputText(key=f"Value_{compound}", visible=True,size=(5,1),
                            default_text=str(fixed_orders_start[compound]))])
                else:
                    layout_fixed_order.append([
                        sg.Text(compound, size=(label_width, 1)),
                        sg.Checkbox("", key=f"Check_{compound}", enable_events=True),
                        sg.InputText(key=f"Value_{compound}", visible=False,size=(5,1),enable_events=True,
                        )])
            # Define the main layout of the calculation settings window. 
            layout_settings = [
            [sg.Text("Calculation Settings", font=("helvetica", 12, "bold"))],
            [sg.Frame("Fit Settings", layout=[
                [sg.Text("Constraint:"),
                sg.Combo(['Monotonic', 'None', 'Via origin'],
                        default_value=calculation_settings["constraint"],
                        key="constraint", enable_events=True, readonly=True,size=(8,1)),
                sg.Text("Degree:"),
                sg.Combo([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], default_value=calculation_settings["deg"], key="deg",size=(3,1))],
                [sg.Text("Goodness-of-fit metric:"),
                sg.Combo(['SE', 'RMSE', 'R2', 'Variance'], default_value=calculation_settings["fit_metric"],
                        key="fit_metric", enable_events=True, readonly=True,size=(8,1))],
            ])],
            [sg.Frame("Order Exploration Settings", layout=[
                [sg.Text("Iterations:"),
                sg.Combo([1, 2, 3, 4, 5, 6, 7,8,9,10,11,12], default_value=calculation_settings["iterations"], key="iterations"),
                sg.Text("Order grid density:"),
                sg.InputText(key="resolution", size=(3, 1),
                            default_text=str(calculation_settings["resolution"]))],
                [sg.Text("Higher initial order mesh density:"),
                sg.Combo(['True','False'], default_value=str(calculation_settings["initial_mesh_denser"]),
                            key='initial_mesh_denser', enable_events=True, readonly=True,size=(6,1))]
            ])],
            [sg.Button("Standard settings",key='standard'),sg.Button("Quick settings",key='quick')],
            [sg.Frame("Order Range to Explore", layout=[
                [sg.Text("From"),
                sg.InputText(key="order_range_low", size=(4, 1),
                            default_text=str(calculation_settings["order_range"][0])),
                sg.Text("to"),
                sg.InputText(key="order_range_high", size=(4, 1),
                            default_text=str(calculation_settings["order_range"][1]))]
            ]),
            sg.Frame("Score Interval", layout=[
                [sg.Text("Interval:"),
                sg.InputText(key="score_interval", size=(4, 1),
                            default_text=str(calculation_settings["score_interval"]))]
            ])],
            [sg.Frame("Fixed Order Species", layout=layout_fixed_order)],
            [sg.Button("OK"), sg.Button("Cancel"),sg.Button("Reset settings")]]

            window_settings = sg.Window("Calculation Settings", layout_settings)
            # Maintain an infinite loop to process events until the window is shut down.
            while True:
                event_settings, values_settings = window_settings.read()
                if event_settings == sg.WIN_CLOSED or event_settings == "Cancel":
                    break
                # Reset the calculation settings if the reset button is pressed. 
                elif event_settings=='Reset settings':
                    # Loop through the elements of the default overlay plot settings and apply these
                    # to the overlay settings window.
                    for i,j in default_calculation_settings.items():
                        if i not in ['order_range','fixed_order_species']:
                            window_settings[i].update(value=str(j))
                    window_settings["order_range_low"].update(value=default_calculation_settings['order_range'][0])
                    window_settings["order_range_high"].update(value=default_calculation_settings['order_range'][1])
                    for compound in columns:
                        window_settings[f"Check_{compound}"].update(value=False)
                        window_settings[f"Value_{compound}"].update(value='')
                        window_settings[f"Value_{compound}"].update(visible=False)
                # Apply the quick settings if the Quick button is pressed. 
                elif event_settings=='quick':
                    window_settings['constraint'].update(value='None')
                    window_settings['resolution'].update(value=4)
                    window_settings['iterations'].update(value=10)
                    window_settings["initial_mesh_denser"].update(value='False')
                # Apply the standard calculation settings if the Standard settings button is pressed.
                elif event_settings=='standard':
                    window_settings['constraint'].update(value='Monotonic')
                    window_settings['resolution'].update(value=7)
                    window_settings['iterations'].update(value=7)
                    window_settings["initial_mesh_denser"].update(value='True')
                # Update the fixed order species frame if the user has ticked or unticked one of the boxes.
                elif event_settings.startswith("Check_"):
                    window_settings[event_settings.replace('Check_','Value_')].update(visible=values_settings[event_settings])
                    if values_settings[event_settings]==False:
                        window_settings[event_settings.replace('Check_','Value_')].update(value='')
                elif event_settings == "OK":
                    # Validate input values.
                    try:
                        for key,value in values_settings.items():
                            if 'Value_' in key:
                                if len(value)>0:
                                    fixed_orders[key.replace('Value_','')]=float(value)
                        order_range_low = float(values_settings["order_range_low"])
                        order_range_high = float(values_settings["order_range_high"])
                        iterations = int(values_settings["iterations"])
                        score_interval = float(values_settings["score_interval"])
                        resolution = float(values_settings["resolution"])
                        deg = int(values_settings["deg"])
                        fit_metric = values_settings["fit_metric"]

                        # Check order range
                        if not (-3 <= order_range_low < order_range_high <= 3.5):
                            raise ValueError("Invalid order range. Must be between -3 and 3.5 with the first value lower than the second.")
                        # Check if the selected fit metric is valid
                        if fit_metric not in ['SE', 'RMSE', 'R2', 'Variance']:
                            raise ValueError("Invalid fit metric. Must be one of: 'SE', 'RMSE', 'R2', 'Variance'.")
                        # Check constraint
                        constraint = values_settings["constraint"]
                        if constraint not in ['Monotonic', 'None', 'Via origin']:
                            raise ValueError("Invalid constraint. Must be 'Monotonic', 'None', or 'Via origin'.")

                        # Check score interval
                        if not (0 <= score_interval <= 100):
                            raise ValueError("Invalid score interval. Must be between 0 and 100.")
                        # Check resolution
                        if not (4 <= resolution <= 30):
                            raise ValueError("Invalid score interval. Must be between 4 and 30.")
                        # Check degree
                        if deg not in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10,11, 12, 13]:
                            raise ValueError("Invalid degree. Must be an integer between 1 and 13.")
                        if constraint=='Via origin':
                            if int(deg)!=1:
                                raise ValueError("Set fitting degree to 1 for linear line via origin.")
                    except ValueError as e:
                        sg.popup_error(f"Invalid input: {e}")
                        continue
                    # Save the calculation settings
                    calculation_settings["constraint"] = constraint
                    calculation_settings["order_range"] = [order_range_low, order_range_high]
                    calculation_settings["iterations"] = round(iterations)
                    calculation_settings["score_interval"] = score_interval
                    calculation_settings["deg"] = round(deg)
                    calculation_settings["fit_metric"] = values_settings["fit_metric"]
                    calculation_settings['initial_mesh_denser']=values_settings['initial_mesh_denser']
                    calculation_settings["resolution"]= round(resolution)
                    calculation_settings["fixed_order_species"]=fixed_orders
                    # Set the fixed order species to None if none have been selected by the user. 
                    if len(fixed_orders)==0:
                        calculation_settings["fixed_order_species"]=None
                    break
            window_settings.close()

        elif event == "Select CSVs":
            # Define counter to know which browse and input text elements to make visible. 
            count=1
            #sg.theme("DefaultNoMoreNagging")
            # Define layout with one set of visible InputText and FilesBrowse elements and several hidden ones. 
            layout = [
                [sg.Text("NB: If >1 file in same browse, alternative names can't be given.\nIf a file is uploaded more than once, the duplicate(s) will be removed.")],
                [sg.InputText(key="file_name_1", enable_events=True, size=(18, 1)),sg.FilesBrowse(key="browse_1"),sg.Text("Data name (optional):",key="text_1"),sg.InputText(key="nickname_1", size=(10, 1))],
                [sg.InputText(key="file_name_2", enable_events=True, size=(18, 1),visible=False), sg.FilesBrowse(key="browse_2",visible=False),sg.Text("Data name (optional):",visible=False,key="text_2"),sg.InputText(key="nickname_2", size=(10, 1),visible=False)],
                [sg.InputText(key="file_name_3", enable_events=True, size=(18, 1),visible=False), sg.FilesBrowse(key="browse_3",visible=False),sg.Text("Data name (optional):",visible=False,key="text_3"),sg.InputText(key="nickname_3", size=(10, 1),visible=False)],
                [sg.InputText(key="file_name_4", enable_events=True, size=(18, 1),visible=False), sg.FilesBrowse(key="browse_4",visible=False),sg.Text("Data name (optional):",visible=False,key="text_4"),sg.InputText(key="nickname_4", size=(10, 1),visible=False)],
                [sg.InputText(key="file_name_5", enable_events=True, size=(18, 1),visible=False), sg.FilesBrowse(key="browse_5",visible=False),sg.Text("Data name (optional):",visible=False,key="text_5"),sg.InputText(key="nickname_5", size=(10, 1),visible=False)],
                [sg.InputText(key="file_name_6", enable_events=True, size=(18, 1),visible=False), sg.FilesBrowse(key="browse_6",visible=False),sg.Text("Data name (optional):",visible=False,key="text_6"),sg.InputText(key="nickname_6", size=(10, 1),visible=False)],
                [sg.InputText(key="file_name_7", enable_events=True, size=(18, 1),visible=False), sg.FilesBrowse(key="browse_7",visible=False),sg.Text("Data name (optional):",visible=False,key="text_7"),sg.InputText(key="nickname_7", size=(10, 1),visible=False)],
                [sg.InputText(key="file_name_8", enable_events=True, size=(18, 1),visible=False), sg.FilesBrowse(key="browse_8",visible=False),sg.Text("Data name (optional):",visible=False,key="text_8"),sg.InputText(key="nickname_8", size=(10, 1),visible=False)],
                [sg.InputText(key="file_name_9", enable_events=True, size=(18, 1),visible=False), sg.FilesBrowse(key="browse_9",visible=False),sg.Text("Data name (optional):",visible=False,key="text_9"),sg.InputText(key="nickname_9", size=(10, 1),visible=False)],
                [sg.InputText(key="file_name_10", enable_events=True, size=(18, 1),visible=False), sg.FilesBrowse(key="browse_10",visible=False),sg.Text("Data name (optional):",visible=False,key="text_10"),sg.InputText(key="nickname_10", size=(10, 1),visible=False)],
                [sg.InputText(key="file_name_11", enable_events=True, size=(18, 1),visible=False), sg.FilesBrowse(key="browse_11",visible=False),sg.Text("Data name (optional):",visible=False,key="text_11"),sg.InputText(key="nickname_11", size=(10, 1),visible=False)],
                [sg.InputText(key="file_name_12", enable_events=True, size=(18, 1),visible=False), sg.FilesBrowse(key="browse_12",visible=False),sg.Text("Data name (optional):",visible=False,key="text_12"),sg.InputText(key="nickname_12", size=(10, 1),visible=False)],
                [sg.InputText(key="file_name_13", enable_events=True, size=(18, 1),visible=False), sg.FilesBrowse(key="browse_13",visible=False),sg.Text("Data name (optional):",visible=False,key="text_13"),sg.InputText(key="nickname_13", size=(10, 1),visible=False)],
                [sg.Button("OK")]
            ]
            # Create the window.
            window_csv = sg.Window("CSV Loader", layout)
            pattern = re.compile(r'/([^/]+)\.csv$')
            # Maintain an infinite loop to process events until the window is shut down.
            while True:
                event, values = window_csv.read()
                if event == sg.WIN_CLOSED:
                    break
                elif event == "OK":
                    try:
                        values_not_nickname=[i[1] for i in list(values.items()) if 'nickname' not in i[0] and i[1]!='']
                        if len(values_not_nickname)==0:
                            break
                        file_names=[]
                        file_paths_full=[]
                        for key,value in values.items():
                            # Loop through all the keys and values in the values of the window_csv.
                            if value==None:
                                break
                            elif "file_name" in key and len(value)!=0:
                                file_paths = value.split(';')
                                file_paths_full=file_paths_full+file_paths
                                # Obtain the index number of the file_name text.
                                index=''.join(char for char in key if char.isdigit())
                                # Obtain the corresponding nickname from the values dictionary. 
                                name=values[f"nickname_{int(index)}"]
                                for i in file_paths:
                                    # Identify the file name from the file path if no nickname has been specified.
                                    if len(name)==0 or len(file_paths)!=1:
                                        name_search = pattern.search(i)
                                        name=name_search.group(1)
                                        file_names.append(name)
                                    else:
                                        file_names.append(name)
                        # Check the file paths and file names lists for duplicates. 
                        # Makes new lists that are guaranteed not to have duplicates. 
                        file_paths_no_duplicates,file_names_no_duplicates=[],[]
                        for count,(i,j) in enumerate(zip(file_paths_full,file_names)):
                            if i not in file_paths_full[:count]:
                                file_paths_no_duplicates.append(i)
                                file_names_no_duplicates.append(j)
                        # Loop through the file paths and create the kinetic data dictionary. 
                        kinetic_data={}
                        for (i,j) in zip(file_paths_no_duplicates,file_names_no_duplicates):
                            dataset=pd.read_csv(f"{i}")
                            # Convert column titles to strings. 
                            dataset.columns = dataset.columns.astype(str)
                            kinetic_data[j]=dataset
                        # Define data as a copy of the kinetic data which is not altered later in the code. 
                        # This allows the original data to be accessed even if the copy is altered. 
                        data = kinetic_data.copy()
                        omission_range_dictionary={}
                        omission_range_dictionary["range_type"]='Percentage'
                        # Figure out if default data cropping mode should be datapoints or range.
                        df_nrows=[]
                        for i in kinetic_data.keys():
                            df_nrows.append(kinetic_data[i].shape[0])
                        if max(df_nrows)>25:
                            omission_mode='Range mode'
                        else:
                            omission_mode='Datapoint mode'
                        # Update the time label and the experiments, column and fixed order species. 
                        time_label=list(data.values())[0].columns[0]
                        experiments = list(data.keys())
                        columns = data[experiments[0]].columns[1:].tolist()  # Exclude the first column
                        # Create a list of reaction species that have concentration profiles that change over time. 
                        reaction_species_not_constant=data[experiments[0]].columns[1:].tolist()
                        for exp in data.keys():
                            for species in columns:
                                # Avoid trying to remove a constant species twice. 
                                if species in reaction_species_not_constant:
                                    # Define the list of concentration values for the relevant species in the relevant experiment. 
                                    values=list(data[exp][species].to_numpy())
                                    # If all values are the same, remove this species from the not constant list
                                    if all([species==values[0] for species in values]):
                                        reaction_species_not_constant.remove(species)
                        fixed_order_species = data[experiments[0]].columns[1:].tolist()  # All reaction species
                        # Update the sheet, columns and result column listboxes with the relevant data
                        # (Sheet names, reactant names and reactant names respectively)
                        window["experiments"].update(values=experiments)
                        window["columns"].update(values=columns)
                        window["result_column"].update(values=reaction_species_not_constant)
                        # Enable various buttons now that the kinetic data has been loaded. 
                        window["Run Auto VTNA"].update(disabled=False)  
                        window["calculation_settings"].update(disabled=False)  
                        window["Check Kinetic Data"].update(disabled=False) 
                        window["Inspect Kinetic Data"].update(disabled=False)
                        window["Generate Overlay Plot"].update(disabled=False)
                        window["Overlay Plot Settings"].update(disabled=False)
                        window["Crop Kinetic Data"].update(disabled=False)
                        files_calculation=file_paths_no_duplicates
                        break
                    except ValueError as e:
                        sg.popup_error(f"Error: {e}")
                        continue
                # If a file is browsed, make the next InputText and FilesBrowser visible. 
                if 'file_name' in event and event[-1]==str(count):
                    file_paths = values[event].split(';')
                    if len(file_paths)>1:
                        window_csv[f"nickname_{count}"].update(visible=False)
                        window_csv[f"text_{count}"].update(visible=False)
                    count+=1
                    window_csv[f"file_name_{count}"].update(visible=True)
                    window_csv[f"browse_{count}"].update(visible=True)
                    window_csv[f"text_{count}"].update(visible=True)
                    window_csv[f"nickname_{count}"].update(visible=True)
                # If the window is closed, collect the filenames selected as a list.
            window_csv.close()

        elif event == "Check Kinetic Data":
            # Apply the check_kinetic_data to investigate whether it has been uploaded correctly 
            # and is ready for analysis. 
            check_kinetic_data(data)

        elif event == "Crop Kinetic Data":
            # Open the kinetic data cropping window either in range mode (for high datapoint
            # density), or in datapoint mode depending the omission_mode variable. 
            if omission_mode=='Range mode':
                data_cropping = kinetic_data.copy()
                # Define a title for each experiment range row. 
                layout_frame_selective=[]
                # Create rows for inputting range crop settings for every experiment data. 
                for key in kinetic_data.keys():
                    row_title = [sg.Text(f'Dataset {key}:',font=("helvetica", 10, "bold"))]
                    # Input boxes for removing data points by applying ranges. 
                    # Populate the InputText boxed with previous range values.
                    if len(omission_range_dictionary)>1:
                        # Create the rows with default values from the previous settings if they exist.
                        remove_data_layout = [
                            sg.InputText(key=f'{key}_START1', size=(5, 1),
                            default_text=omission_range_dictionary[f'{key}_START1']),
                            sg.Text('to time'),
                            sg.InputText(key=f'{key}_END__1', size=(5, 1),
                            default_text=omission_range_dictionary[f'{key}_END__1']),
                            sg.Text('and from time'),
                            sg.InputText(key=f'{key}_START2', size=(5, 1),
                            default_text=omission_range_dictionary[f'{key}_START2']),
                            sg.Text('to time'),
                            sg.InputText(key=f'{key}_END__2', size=(5, 1),
                            default_text=omission_range_dictionary[f'{key}_END__2']),
                            sg.Text("Keep only every:"),
                            sg.InputText(key=f'{key}_less_dense', enable_events=True,size=(4,1),
                            default_text=omission_range_dictionary[f'{key}_less_dense']),
                            sg.Text("datapoint.")]
                    else:
                        # Input boxes for removing data points by applying ranges (no values 
                        # already saved)
                        remove_data_layout = [
                            sg.Text(f'Remove datapoints from time '),
                            sg.InputText(key=f'{key}_START1', size=(5, 1)),
                            sg.Text('to time'),
                            sg.InputText(key=f'{key}_END__1', size=(5, 1)),
                            sg.Text('and from time'),
                            sg.InputText(key=f'{key}_START2', size=(5, 1)),
                            sg.Text('to time'),
                            sg.InputText(key=f'{key}_END__2', size=(5, 1)),
                            sg.Text("Keep only every:"),
                            sg.InputText(key=f'{key}_less_dense', enable_events=True,size=(4,1)),
                            sg.Text("datapoint.")]
                    RM_line=row_title+remove_data_layout
                    layout_frame_selective.append(RM_line)
                # Create the frame for applying cropping to each experiment separately.
                remove_data_selective=[sg.Frame("Apply cropping separately for each experiment data:",layout=
                layout_frame_selective,visible=True,key='selective_frame')]
                # Create a list containing a frame for applying cropping to all experiments at once. 
                if len(omission_range_dictionary)>1:
                    # Create the frame with default values previously selected. 
                    remove_data_total=[sg.Frame("Apply cropping to all experiment data at once:",layout=[
                            [sg.Text("NB: Only applied to datasets for which settings haven't \
    been specified above.")],[sg.Text(f'Remove datapoints from time '),
                            sg.InputText(key='total_START1', size=(5, 1),
                            default_text=omission_range_dictionary['total_START1']),
                            sg.Text('to time'),
                            sg.InputText(key='total_END__1', size=(5, 1),
                            default_text=omission_range_dictionary['total_END__1']),
                            sg.Text('and from time'),
                            sg.InputText(key='total_START2', size=(5, 1),
                            default_text=omission_range_dictionary['total_START2']),
                            sg.Text('to time'),
                            sg.InputText(key='total_END__2', size=(5, 1),
                            default_text=omission_range_dictionary['total_END__2']),
                            sg.Text("Keep only every:"),
                            sg.InputText(key='total_less_dense', size=(4,1),
                            default_text=omission_range_dictionary['total_less_dense']),
                            sg.Text("datapoint.")
                            ]],visible=True,key='total_frame')]
                
                else:
                    # Create the frame with no default values as none have previously been selected. 
                    remove_data_total=[sg.Frame("Apply cropping to all experiment data at once:",layout=[
                            [sg.Text("NB: Only applied to datasets for which settings haven't \
    been specified above.")],[sg.Text(f'Remove datapoints from time '),
                            sg.InputText(key='total_START1', size=(5, 1)),
                            sg.Text('to time'),
                            sg.InputText(key=f'total_END__1', size=(5, 1)),
                            sg.Text('and from time'),
                            sg.InputText(key='total_START2', size=(5, 1)),
                            sg.Text('to time'),
                            sg.InputText(key='total_END__2', size=(5, 1)),
                            sg.Text("Keep only every:"),
                            sg.InputText(key='total_less_dense', enable_events=True,size=(4,1)),
                            sg.Text("datapoint.")
                            ]],visible=True,key='total_frame')]
                # Complete the layout by combining the experiment specific and total layout lists. 
                rows_range=[remove_data_selective]+[remove_data_total]
                # Add the final elements to the range cropping layout.
                rows_range.append([sg.Text("Relative or absolute times"),
                                    sg.Combo(['Absolute','Percentage'], default_value=omission_range_dictionary["range_type"],
                                            key="range_type", enable_events=True, readonly=True)])
                layout_range = rows_range + [[sg.Button('OK'),sg.Button("Reset")]]
                # Create the range crop window.
                window_range = sg.Window('Crop Data Range', layout_range)
                # Maintain an infinite loop to process events until the window is shut down.
                while True:
                    event_range, values_range = window_range.read()
                    if event_range=="Reset":
                        data_cropping=kinetic_data.copy()
                        for key in kinetic_data.keys():
                            window_range[f'{key}_START1'].update(value='')
                            window_range[f'{key}_START2'].update(value='')
                            window_range[f'{key}_END__1'].update(value='')
                            window_range[f'{key}_END__2'].update(value='')
                        window_range['total_START1'].update(value='')
                        window_range['total_START2'].update(value='')
                        window_range['total_END__1'].update(value='')
                        window_range['total_END__2'].update(value='')
                    if event_range == 'OK':
                        try:
                            # Create a values dictionary without empty values.
                            values_range_cleaned={}
                            for key,value in values_range.items():
                                if value!=None and value!='' and key!='range_type':
                                    values_range_cleaned[key]=value
                            range_info=extract_range_info(values_range_cleaned,data)
                            # Check if all inputs are numerical:
                            num_check=[is_float(i) for i in list(values_range_cleaned.values())]
                            if not all(num_check):
                                raise ValueError("Not all the range values provided are numerical.")
                            # Check if range value inputs are correct.
                            if check_ends(values_range_cleaned,data)==False:
                                raise ValueError("Invalid range inputs. Make sure that each connected \
    start and end value is provided for each row.")
                            if not check_increase(range_info):
                                raise ValueError("Some of the start and end range value pairs are not increasing.")
                            for i in range_info.keys():
                                if len(range_info[i])==4:
                                    if float(range_info[i][1])==100:
                                        raise ValueError(f"The first range interval for {i} can't reach to 100. \
    Make this the second interval rather.")
                        except ValueError as e:
                            sg.popup_error(f"Invalid input: {e}")
                            continue
                        if True in ['dense' in i for i in list(values_range_cleaned.keys())]:
                            # Remove datapoints from data_cropping at regular intervals inputed:
                            selection_dictionary_lower_density=auto_vtna.make_VTNA_selection(kinetic_data).copy()
                            sparser_total=1
                            # If the total less dense window has been populated by the user, define the total sparser factor. 
                            if 'total_less_dense' in list(values_range_cleaned.keys()):
                                sparser_total=abs(int(float(values_range_cleaned['total_less_dense'])))
                            # Apply data sparsing to each experiment.
                            for i in selection_dictionary_lower_density['RM'].keys():
                                sparser=sparser_total
                                if f'{i}_less_dense' in values_range_cleaned.keys():
                                    sparser=abs(int(float(values_range_cleaned[f'{i}_less_dense'])))
                                length=len(kinetic_data[i].iloc[:, 0])
                                all_DP=list(range(0,length))
                                # Define a list of datapoints to keep.
                                to_keep=all_DP[0::sparser]
                                # If sparsing removed the last value, add it back. 
                                if to_keep[-1]!=all_DP[-1]:
                                    to_keep.append(all_DP[-1])
                                # Define which datapoints need to be removed.
                                to_remove=list(set(all_DP) - set(to_keep))
                                # Don't define omissions datapoints if none are to be removed.
                                if len(to_remove)==0:
                                    continue
                                # Update the VTNA selection dictionary with the list of datapoints to be removed. 
                                selection_dictionary_lower_density['RM'][i]['omissions']=to_remove
                            # Apply the VTNA selection dictionary containing lists of datapoints to be removed for each 
                            # dataset using the VTNA_omissions function. 
                            data_cropping=VTNA_omissions(data_cropping,selection_dictionary_lower_density)
                        # make all values in range_info numerical.
                        range_info={key: [float(value) for value in values] for key, values in range_info.items()}
                        # Prepare a new VTNA selection dictionary to apply range omissions.
                        selection_dictionary_ranges=auto_vtna.make_VTNA_selection(kinetic_data).copy()
                        range_list=['range']
                        # Save the range omissions settings for each experiment
                        for i in kinetic_data.keys():
                            if i in list(range_info.keys()):
                                selection_dictionary_ranges['RM'][i]['omissions']=range_list+range_info[i][:2]
                            else:
                                if 'total' in list(range_info.keys()):
                                    selection_dictionary_ranges['RM'][i]['omissions']=range_list+range_info['total'][:2]
                        # Apply the range omissions by calling VTNA_omissions. 
                        data_cropping=VTNA_omissions(data_cropping,selection_dictionary_ranges,range_mode=values_range["range_type"])
                        # Check if there are any rows in the range cropping window with more than one range defined. 
                        if 4 in [len(i) for i in range_info.values()]:
                            selection_dictionary_ranges_2=auto_vtna.make_VTNA_selection(kinetic_data).copy()
                            # Loop through the experiments i to apply the second omission range if selected by the user. 
                            for i in kinetic_data.keys():
                                if i in list(range_info.keys()):
                                    if len(range_info[i])==4:
                                        selection_dictionary_ranges_2['RM'][i]['omissions']=range_list+range_info[i][2:]
                                # If no second omissions range has been defined for a given dataset, apply the total 
                                # second omission range if present. 
                                else:
                                    if 'total' in list(range_info.keys()):
                                        if len(range_info['total'])==4:
                                            selection_dictionary_ranges_2['RM'][i]['omissions']=range_list+range_info['total'][2:]
                            # Apply the second range omissions by calling VTNA_omissions.
                            data_cropping=VTNA_omissions(data_cropping,selection_dictionary_ranges_2)
                        # Update the data variable to equal the cropped version. 
                        data=data_cropping
                        # Save the range values so that they can be shown next time a range cropping window is generated. 
                        omission_range_dictionary=values_range
                        window_range.close()
                    if event_range == sg.WINDOW_CLOSED:
                        window_range.close()
                        break
            # Apply datapoint mode data cropping. 
            else:
                rows = []
                # Create the layout for the datapoint cropping window using the previously saved information. 
                if omission_dictionary!=None:
                    for key, df in kinetic_data.items():
                        row_title = [sg.Text(f'Row {key}')]
                        checkboxes = [sg.Checkbox(str(i), key=f'CHECK_{key}__{i}',default=omission_dictionary[f'CHECK_{key}__{i}']) for i in range(1, len(df[time_label]) + 1)]
                        rows.append(row_title + checkboxes)
                else:
                # Create the layout for the datapoint cropping window if no information has already been saved.
                    for key, df in kinetic_data.items():
                        row_title = [sg.Text(f'Row {key}')]
                        checkboxes = [sg.Checkbox(str(i), key=f'CHECK_{key}__{i}') for i in range(1, len(df[time_label]) + 1)]
                        rows.append(row_title + checkboxes)
                layout = rows + [[sg.Button('OK'),sg.Button("Reset")]]
                window_omissions = sg.Window('Select rows to reversibly remove from kinetic data.', layout)
                # Maintain an infinite loop to process events until the window is shut down.
                while True:
                    event_omissions, values_omissions = window_omissions.read()
                    # Remove all box ticks if the reset button is pressed by the user. 
                    if event_omissions=="Reset":
                        for key, df in kinetic_data.items():
                            for i in range(1, len(df[time_label])):
                                window_omissions[f'CHECK_{key}__{i}'].update(value=False)
                    omission_dictionary=values_omissions
                    selection_dictionary=auto_vtna.make_VTNA_selection(kinetic_data)
                    # Create a copy of the selection dictionary and values to the relevant subdictionary 
                    # corresponding to the boxes that have been checked. 
                    dic=selection_dictionary.copy()
                    for i in dic['RM'].keys():
                        rm_list=[]
                        for j in values_omissions.keys():
                            if i in j and values_omissions[j]:
                                rm_list.append(int(j[-3:].replace('_',''))-1)
                        if len(rm_list)>0:
                            dic['RM'][i]['omissions']=rm_list
                    # Apply the omissions to the kinetic data using the "dic" VTNA selection dictionary. 
                    data=VTNA_omissions(kinetic_data,dic)
                    if event_omissions == sg.WINDOW_CLOSED or event_omissions == 'OK':
                        window_omissions.close()
                        break

        elif event == "Inspect Kinetic Data":
            # Open a new window for inspecting kinetic data
            # Ensure that all exisiting plots are closed to avoid bugs.
            plt.close('all')
            # Define the layout of the window including both buttons for different plots and 
            # the plot settings. 
            layout_inspect_kinetic_data = [
                [sg.Button("Generate initial concentration table")],
                [sg.Button("Plot kinetic data")],
                [sg.Button("Plot data for selected experiments & species")],
                [sg.Button("Plot concentration profiles")],
                [sg.Button('Export the kinetic dataset to .xlsx or .csv')],
                [sg.Frame("Data Visualisation Settings", layout=[
                    [sg.Text("Concentration unit:"),
                     sg.InputText(key="y_unit", size=(5, 1),default_text=str(data_plotting_settings["y_unit"]),enable_events=True),
                     sg.Text("Time unit: "),
                     sg.InputText(key="t_unit", size=(5, 1),default_text=str(data_plotting_settings["t_unit"]),enable_events=True)],
                    [sg.Text("Figure size scaler: "),
                     sg.InputText(key="size_scaler", size=(5, 1),default_text=float(data_plotting_settings["size_scaler"]),enable_events=True),
                     sg.Text("Max y:     "),
                     sg.InputText(key="ylim", size=(5, 1),default_text=str(data_plotting_settings["ylim"]),enable_events=True)],
                    [sg.Text("Datapoint scaler:   "),
                     sg.InputText(key="DP_scaler", size=(5, 1),default_text=float(data_plotting_settings["DP_scaler"]),enable_events=True),
                     sg.Text("Linewidth:"),
                     sg.InputText(key="linewidth", size=(5, 1),default_text=float(data_plotting_settings["linewidth"]),enable_events=True)],
                    [sg.Text("Legend position:    "),
                     sg.Combo(['Inside','Outside'], size=(7, 1),default_value=data_plotting_settings["legend_position"],key="legend_position", enable_events=True, readonly=True)],
                    [sg.Text("Significant figures for initial conc. table:"),
                     sg.InputText(key="SF", size=(5, 1),default_text=float(data_plotting_settings["SF"]),enable_events=True)],
                    [sg.Text("Mode for visualising kinetic data:"),sg.Combo(['Scroll', 'Together', 'Separate'], default_value=data_plotting_settings["mode"],
                            key="mode", enable_events=True, readonly=True)],
                ])],
                [sg.Button("Back"),sg.Button('Reset settings')]
            ]
            # Create the window. 
            window_inspect_kinetic_data = sg.Window("Inspect Kinetic Data", layout_inspect_kinetic_data)
            # Maintain an infinite loop to process events until the window is shut down.
            while True:
                event_inspect_kinetic_data, values_inspect_kinetic_data = window_inspect_kinetic_data.read()
                if event_inspect_kinetic_data == sg.WIN_CLOSED or event_inspect_kinetic_data == "Back":
                    break
                # Reset the settings to the default values if the user presses reset settings. 
                elif event_inspect_kinetic_data=='Reset settings':
                    # Loop through the elements of the default data plotting settings and apply these
                    # to the inspect kinetic data window.
                    for i,j in default_data_plotting_settings.items():
                        if 'scaler' in i or i in ['linewidth','SF']:
                            j=float(j)
                        else:
                            j=str(j)
                        window_inspect_kinetic_data[i].update(value=j)
                # Check that the inputed settings are valid before creating any plots or table. 
                if event_inspect_kinetic_data in ["Generate initial concentration table","Plot kinetic data",
                        "Plot data for selected experiments & species","Plot concentration profiles"]:
                    try:
                        if not is_float(values_inspect_kinetic_data["ylim"]):
                            ylim = None
                        else: 
                            ylim = float(values_inspect_kinetic_data["ylim"])
                        y_unit = values_inspect_kinetic_data["y_unit"]
                        t_unit = values_inspect_kinetic_data["t_unit"]
                        if not is_float(values_inspect_kinetic_data["SF"]):
                            raise ValueError("Significant figures for initial concentration table must be a numerical value.")
                        significant_figures=int(abs(float(values_inspect_kinetic_data["SF"])))
                        data_plotting_settings["legend_position"]=values_inspect_kinetic_data["legend_position"]
                        DP_scaler = float(values_inspect_kinetic_data["DP_scaler"])
                        linewidth = float(values_inspect_kinetic_data["linewidth"])
                        mode = values_inspect_kinetic_data["mode"]
                        if len(t_unit)==0 or t_unit=='None':
                            t_unit=None
                        size_scaler = float(values_inspect_kinetic_data["size_scaler"])
                        # Check that inputs are correctly defined to avoid errors. 
                        # Generate error pop-up windows if there are issues. 
                        if not (0.1 < size_scaler <= 3):
                            raise ValueError("Invalid figure size scaler. Must be between 0.1 and 3.")
                        if not (0 <= DP_scaler <= 5):
                            raise ValueError("Invalid datapoint size scaler. Must be between 0 and 5.")
                        if not (0 <= linewidth <= 5):
                            raise ValueError("Invalid linewidth. Must be between 0 and 5.")
                        if type(ylim)==float:
                            if not (0 < ylim):
                                raise ValueError("Invalid concentration axis limit. Must be a number above 0.")
                        if type(ylim)==str:
                            if ylim!='None' and ylim!='':
                                raise ValueError("Concentration axis limit must be a numerical value.")
                        if significant_figures==0 or significant_figures>9:
                            raise ValueError("Significant figures for initial concentration must be between 1 and 9.")
                        if type(ylim)==str:
                            if ylim!='None' and ylim!='':
                                raise ValueError("Concentration axis limit must be a numerical value.")
                        data_plotting_settings["SF"]=int(significant_figures)
                        data_plotting_settings["size_scaler"]=size_scaler
                        # Ensure that the ylim setting is set to None rather than '' or 'None' if this is selected. 
                        if ylim not in ['None','']:
                            data_plotting_settings["ylim"]=ylim
                        else:
                            data_plotting_settings["ylim"]=None
                        data_plotting_settings["y_unit"]=y_unit
                        data_plotting_settings["t_unit"]=t_unit
                        data_plotting_settings["mode"]=mode
                        data_plotting_settings["DP_scaler"]=DP_scaler
                        data_plotting_settings["linewidth"]=linewidth
                    except ValueError as e:
                            sg.popup_error(f"Invalid input: {e}")
                            continue
                # React to the specifics of the users button click. 
                if event_inspect_kinetic_data == "Generate initial concentration table":
                    generate_initial_concentration_table(data,size_scaler,y_unit,significant_figures=data_plotting_settings["SF"])
                elif event_inspect_kinetic_data == "Plot kinetic data":
                    # Ask the user whether to visualise mass balance using a custom popup window.
                    layout_visualise_mass_balance = [
                        [sg.Text("Do you want to visualise mass balance for the kinetic data?")],
                        [sg.Button("Yes"), sg.Button("No")]
                    ]
                    window_visualise_mass_balance = sg.Window("Visualise Mass Balance", layout_visualise_mass_balance)
                    event_visualise_mass_balance, values_visualise_mass_balance = window_visualise_mass_balance.read()
                    window_visualise_mass_balance.close()
                    if event_visualise_mass_balance==sg.WIN_CLOSED:
                        continue
                    if event_visualise_mass_balance == "Yes":
                        # Use Listbox to create a multi-selection list
                        selected_species = sg.Listbox(
                            values=columns,
                            select_mode=sg.LISTBOX_SELECT_MODE_MULTIPLE,
                            size=(20, min(12, len(fixed_order_species))),
                            key="selected_species_listbox"
                        )
                        # Create a new popup window to enable the user to select the reaction 
                        # species for the mass balance calculations. 
                        layout_select_species = [
                            [sg.Text("Select reaction species for mass balance calculations:")],
                            [selected_species],
                            [sg.Button("OK"), sg.Button("Cancel")]
                        ]
                        window_select_species = sg.Window("Select Reaction Species", layout_select_species)
                        # Maintain an infinite loop to process events until the window is shut down.
                        while True:
                            event_select_species, values_select_species = window_select_species.read()
                            # Close the window if the user closes it or clicks cancel. 
                            if event_select_species == sg.WIN_CLOSED or event_select_species == "Cancel":
                                window_select_species.close()
                                break
                            elif event_select_species == "OK":
                                selected_species_values = values_select_species["selected_species_listbox"]
                                if not selected_species_values:
                                    sg.popup_error("Please select one or more reaction species.")
                                else:
                                    # Query user for stoichiometry numbers.
                                    layout_stoichiometry_values = [
                                        [sg.Text(f"Enter stoichiometry values for selected reaction species:")],
                                    ]
                                    # Define columns with reaction species and with input boxes for order values:
                                    species_column=[]
                                    orders_column=[]
                                    for species in selected_species_values:
                                        species_column.append([sg.Text(f"{species}")])
                                        orders_column.append([sg.InputText(key=f"stoichiometry_value_{species}", size=(5, 1))])
                                    # Add the columns to the layout. 
                                    layout_stoichiometry_values.append([sg.Column(species_column,element_justification='left'),sg.Column(orders_column,element_justification='left')])
                                    layout_stoichiometry_values.append([sg.Button("OK"), sg.Button("Cancel")])
                                    # Create the popup window for inputing stochiometry values. 
                                    window_stoichiometry_values = sg.Window("Enter Stoichiometry Values", layout_stoichiometry_values)
                                    # Maintain an infinite loop to process events until the window is shut down.
                                    while True:
                                        event_stoichiometry_values, values_stoichiometry_values = window_stoichiometry_values.read()
                                        if event_stoichiometry_values == sg.WIN_CLOSED or event_stoichiometry_values == "Cancel":
                                            break
                                        elif event_stoichiometry_values == "OK":
                                            try:
                                                # Save the stochiometry values inputted by the user as a dictionary. 
                                                stoichiometry_values = {}
                                                for species in selected_species_values:
                                                    stoichiometry_value = values_stoichiometry_values.get(f"stoichiometry_value_{species}", "")
                                                    stoichiometry_values[species] = float(stoichiometry_value)
                                                # Call the plot kinetic data functio using the stochiometry values dictionary. 
                                                plot_kinetic_data(data, stoichiometry_values,data_plotting_settings)
                                                break
                                            except ValueError as e:
                                                sg.popup_error(f"Invalid input: {e}")
                                    window_stoichiometry_values.close()
                                break
                        window_select_species.close()
                    else:
                        # If the user chooses not to visualise mass balance, call the plot_kinetic_data function without 
                        plot_kinetic_data(data,None,data_plotting_settings)
                elif event_inspect_kinetic_data == "Plot data for selected experiments & species":  
                    # Selected the experiments and reaction species selected by the user. 
                    selected_experiments = values["experiments"]
                    selected_species = values["columns"]
                    # Create an error message is no selected experiments or species have been selected. 
                    if not selected_experiments and not selected_species:
                        sg.popup_error("No experiments nor reaction species have been selected.")
                    elif not selected_experiments:
                        sg.popup_error("No experiments have been selected.")
                    elif not selected_species:
                        sg.popup_error("No reaction species have been selected.")
                    else:
                        # Create a modified kinetic data variable. First, add the dataframes for each selected sheet.
                        data_2=copy.deepcopy(data)
                        data_modified={}
                        for exp in selected_experiments:
                            data_modified[exp]=data_2[exp]
                        # Remove the concentration profiles not included in "selected_species" from the kinetic data.
                        for RM in data_modified.keys():
                            for i, NS in enumerate(data_modified[RM]):
                                if NS not in selected_species and i!=0:
                                    del data_modified[RM][NS]
                        # Ask the user whether to visualise mass balance using a custom popup window. 
                        layout_visualise_mass_balance = [
                            [sg.Text("Do you want to visualise mass balance for the kinetic data?")],
                            [sg.Button("Yes"), sg.Button("No")]]
                        # Create the visualise mass balance window. 
                        window_visualise_mass_balance = sg.Window("Visualise Mass Balance", layout_visualise_mass_balance)
                        event_visualise_mass_balance, values_visualise_mass_balance = window_visualise_mass_balance.read()
                        window_visualise_mass_balance.close()
                        if event_visualise_mass_balance == sg.WIN_CLOSED:
                            continue
                        if event_visualise_mass_balance == "Yes":
                            # Use Listbox to create a multi-selection list
                            selected_species_list = sg.Listbox(
                                values=columns,
                                select_mode=sg.LISTBOX_SELECT_MODE_MULTIPLE,
                                size=(20, min(12, len(selected_species))),
                                key="selected_species_listbox"
                            )
                            layout_select_species = [
                                [sg.Text("Select reaction species for mass balance calculations:")],
                                [selected_species_list],
                                [sg.Button("OK"), sg.Button("Cancel")]
                            ]
                            window_select_species = sg.Window("Select Reaction Species", layout_select_species)
                            # Maintain an infinite loop to process events until the window is shut down.
                            while True:
                                event_select_species, values_select_species = window_select_species.read()
                                if event_select_species == sg.WIN_CLOSED or event_select_species == "Cancel":
                                    window_select_species.close()
                                    break
                                elif event_select_species == "OK":
                                    selected_species_values = values_select_species["selected_species_listbox"]
                                    if not selected_species_values:
                                        sg.popup_error("Please select one or more reaction species.")
                                    else:
                                        # Query user for stoichiometry numbers.
                                        layout_stoichiometry_values = [[sg.Text(f"Enter stoichiometry values for selected reaction species:")]]
                                       # Define columns with reaction species and with input boxes for order values:
                                        species_column=[]
                                        orders_column=[]
                                        for species in selected_species_values:
                                            species_column.append([sg.Text(f"{species}")])
                                            orders_column.append([sg.InputText(key=f"stoichiometry_value_{species}", size=(5, 1))])
                                        # Add the columns to the layout. 
                                        layout_stoichiometry_values.append([sg.Column(species_column,element_justification='left'),sg.Column(orders_column,element_justification='left')])
                                        layout_stoichiometry_values.append([sg.Button("OK"), sg.Button("Cancel")])
                                        # Create the window for the user to input stochiometry values for mass balance.
                                        window_stoichiometry_values = sg.Window("Enter Stoichiometry Values", layout_stoichiometry_values)
                                        # Maintain an infinite loop to process events until the window is shut down.
                                        while True:
                                            event_stoichiometry_values, values_stoichiometry_values = window_stoichiometry_values.read()
                                            if event_stoichiometry_values == sg.WIN_CLOSED or event_stoichiometry_values == "Cancel":
                                                break
                                            elif event_stoichiometry_values == "OK":
                                                stoichiometry_values = {}
                                                # Save the selected stochiometry values. 
                                                for species in selected_species_values:
                                                    stoichiometry_value = values_stoichiometry_values.get(f"stoichiometry_value_{species}", "")
                                                    stoichiometry_values[species] = stoichiometry_value
                                                # Plot the modified kinetic data with the inputted stochiometry values and the data plotting settings. 
                                                plot_kinetic_data(data_modified, stoichiometry_values,data_plotting_settings)
                                                break
                                        window_stoichiometry_values.close()             
                            window_select_species.close()
                        else:
                            # If the user chooses not to visualise mass balance, return None for mass balance and stoichiometry list
                            plot_kinetic_data(data_modified,None,data_plotting_settings)

                elif event_inspect_kinetic_data == "Plot concentration profiles":
                    # Create a list of reaction species 
                    reaction_species_list = fixed_order_species
                    # Create a list of selected experiments:
                    selected_experiments = values["experiments"]
                    # Use Combo to create a dropdown menu for the user to select reaction species. 
                    selected_species = sg.Combo(
                        values=reaction_species_list,
                        default_value=reaction_species_list[0] if reaction_species_list else "",
                        size=(20, 1),
                        key="selected_species_dropdown",)
                    # Create the select species layout. 
                    layout_select_species = [
                        [sg.Text("Select one reaction species to plot concentration profiles:")],
                        [selected_species],
                        [sg.Button("OK"), sg.Button("Cancel")]]
                    if len(selected_experiments)>0:
                        layout_select_species.insert(2,[sg.Checkbox("Only include selected experiments",key='experiments',default=Plot_concentration_profiles_selected_experiments,enable_events=True)])
                    else:
                        layout_select_species.insert(2,[sg.Text('NB: Select experiments in main window for subset of experiments.', font=('Arial', 9))])
                    # Create a window for selecting reaction species. 
                    window_select_species = sg.Window("Select Reaction Species", layout_select_species)
                    # Maintain an infinite loop to process events until the window is shut down.
                    while True:
                        event_select_species, values_select_species = window_select_species.read()
                        if event_select_species == sg.WIN_CLOSED or event_select_species == "Cancel":
                            break
                        elif event_select_species == "OK":
                            selected_species_value = values_select_species["selected_species_dropdown"]
                            if len(selected_experiments)==0:
                                plot_concentration_profiles(data, selected_species_value,data_plotting_settings)
                            else:
                                if values_select_species['experiments']:
                                    # Make a copy of the data and make a modified version of it with only the selected experiments. 
                                    data_2=copy.deepcopy(data)
                                    data_modified={}
                                    for exp in selected_experiments:
                                        data_modified[exp]=data_2[exp]
                                    selected_species_value = values_select_species["selected_species_dropdown"]
                                    # Plot the modified kinetic data. 
                                    plot_concentration_profiles(data_modified, selected_species_value,data_plotting_settings)
                                else:
                                    plot_concentration_profiles(data, selected_species_value,data_plotting_settings)
                                break
                    window_select_species.close()
                elif event_inspect_kinetic_data=='Export the kinetic dataset to .xlsx or .csv':
                        # Open a browse window where the user can save the active kinetic data ("data").
                        file_path_save_data = sg.popup_get_file('Save As', save_as=True, no_window=True,file_types=(("Excel File", "*.xlsx"), ("CSV File", "*.csv")))
                        # Identify the selected extension (xlsx or csv)
                        file_extension = os.path.splitext(file_path_save_data)[1]
                        if file_extension=='.xlsx':
                            # Save an .xlsx file with one experiment per sheet. 
                            with pd.ExcelWriter(file_path_save_data) as writer:
                                for name, df in data.items():
                                    df.to_excel(writer, sheet_name=name,index=False)
                        elif file_extension=='.csv':
                            # Add a confirmation box to inform the user that saving as csv will generate one file per experiment. 
                            proceed = sg.popup_ok_cancel('NB: Saving as CSV will generate one file per experiment in the selected folder.\nDo you want to proceed?', title='Confirmation')
                            # Save each of the experiments as CSVs if the user has clicked OK. 
                            if proceed == 'OK':
                                for name, df in data.items():
                                    df.to_csv(os.path.join(os.path.dirname(file_path_save_data), f"{name}.csv"), index=False)
                            else:
                                continue
            window_inspect_kinetic_data.close()
    window.close()
if __name__ == "__main__":
    main()

