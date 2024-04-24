from auto_vtna.VTNA_functions import VTNA_orders, on_click
import numpy as np
import pandas as pd
import itertools
from functools import partial
import matplotlib.pyplot as plt
from openpyxl.drawing.image import Image
from io import BytesIO
import matplotlib.ticker as tkr
import warnings

class Automatic_VTNA:
    def __init__(self,data,VTNA_selection,order_range=[-1.5,2.5],resolution=10,deg=5, \
                 fit_metric='SE',iterations=4,constraint='monotonic',score_interval=0.15,\
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
        self.results, self.best_orders, self.interval, self.best_slope = VTNA_orders(self.data,self.VTNA_selection, \
        self.order_range,self.resolution,self.deg,self.fit_metric,self.iterations,self.constraint, \
        self.score_interval,self.fixed_order_species,self.initial_mesh_denser)

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
                    fixed_order_string=fixed_order_string+f"\nFixed order for {i}:{self.VTNA_selection['normalised_species'][i]}."
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

