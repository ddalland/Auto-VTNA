from auto_vtna.VTNA_functions import VTNA_new,data_scaling_y,score_fit,origin_line,simple_line
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from polyfit import PolynomRegressor, Constraints
import warnings

class Normal_VTNA:
    def __init__(self,data,VTNA_selection):
        # Initialise the Normal VTNA class. 
        self.data=data
        self.VTNA_selection=VTNA_selection
        # Perform normal VTNA calculation:
        self.result= VTNA_new(self.data,self.VTNA_selection)
        # Check if the concentration profiles begin at 0 (products) or not (reactants).
        if 'RM' in list(VTNA_selection.keys()):
            initial_output_values=[data[RM][VTNA_selection['output_species']].to_numpy()[0] for RM in list(VTNA_selection['RM'].keys())]
        else: 
            initial_output_values=[data[RM][VTNA_selection['output_species']].to_numpy()[0] for RM in list(data.keys())]
        if all(elem==initial_output_values[0] for elem in initial_output_values):
            self.shifted_output=False
        else:
            self.shifted_output=True
    def plot_VTNA(self,y_unit='M',size_scaler=1,title=None,grid=False, \
                  xaxis_size_scaler=None,DP_scaler=1,xtick_rotation=0,
                  show_legend=True,legend_outside=False,plt_show=True,
                  linewidth=1,scale_xaxis=False,xaxis_notation='Automatic'):
        """
        Plots the VTNA overlay plot for the kinetic data according to the specifications of the \
        VTNA selection dictionary. 
        Args:
            y_unit (str): String defining the output concentration unit for the overlay plots.
            size_scaler (int or float): Factor by which the overlay plot size is adjusted.
            title (str): String defining the title of the plot. 
            grid (bol): Set to True to create a grid based on x and y axis tick locations. 
            xaxis_size_scaler (int or float): Factor by which the size of the xaxis label is adjusted. 
            xtick_rotation (int or float): Determines the degree of rotation of the xaxis tick numbers. 
            show_legend (bol): Set to True to show the legend for different reaction profiles or to False \
            to omit it. 
            legend_outside (bol): Set to True to show the legend for the different reaction profiles outside \
            the main plot. Otherwise, it will be placed automatically at the most open space within the plot \
            automatically. 
            plt_show (bol): Set to True to apply plt.show() at the end of the code to ensure that the plot is \
            generated. Otherwise if set to False, the figure can be modified externally after the execution of the code, giving \
            greater custumisability. 
            linewidth (float): Width of the line that connects the datapoints of concentration profiles. 
            scale_xaxis (bol): Set to True to scale the x axis values (normalised time axis) from 0 to 1. 
            xaxis_notation (string): Determines the notation of the x axis. Can be set to 'Normal' to avoid scientific \
            notation, 'Scientific' to ensure scienfific notation eg. 1e2 instead of 100, or 'Automatic' to get the default setting. 
        """
        plt.rcParams['font.family'] = 'DejaVu Sans'
        # Determine the shape of the overlay plot
        fig_size=(5*size_scaler,4*size_scaler)
        fig=plt.figure(figsize=fig_size)
        # Set the datapoint size scaler to 1 if an incorrect value has been assigned.
        if DP_scaler==None or isinstance(DP_scaler, str):
            DP_scaler=1
        # Determine the normalised and output reaction species from the selection dictionary.
        output_species=self.VTNA_selection['output_species']
        normalised_species=[]
        # Create a list of the order values of the different normalised reaction speices. 
        orders=[]
        for ns in self.VTNA_selection['normalised_species'].keys():
            normalised_species.append(ns)
            orders.append(str(self.VTNA_selection['normalised_species'][ns]))
        # If scale_xaxis is set to true, scale the transformed time values by the largest across the dataset. 
        scaler=1
        if scale_xaxis:
            tTs=[]
            for i in self.result.keys():
                tTs=tTs+list(self.result[i]['tT'])
            scaler=max(tTs)
        # Loop through each experiment in the result from the normal VTNA calculation. 
        for exp in self.result.keys():
            # Plot the overlay curve for each experiment dataset i
            x_values=[float(j)/scaler for j in self.result[exp]['tT'].to_numpy()]
            ys=[float(j) for j in self.result[exp][output_species].to_numpy()]
            plt.plot(x_values,ys,marker='o',label=exp,markersize=6*DP_scaler,linewidth=linewidth)
        if show_legend:
            if legend_outside:
                # Calculate the width of the legend dynamically based on the longest text element
                legend = plt.legend()
                # Obtain the legend width in display units. 
                legend_bbox = legend.get_window_extent()
                legend_width = legend_bbox.width
                # Obtain the figure width in display units.
                fig_box = plt.gcf().get_window_extent()
                fig_width_display = fig_box.width
                # See by what fraction the figure width needs to be extended for outside legend.
                legend_width_scaled=legend_width/fig_width_display
                # Scale the right argument according to the legend width
                plt.subplots_adjust(right=1 - legend_width_scaled)
                # Adjust top and bottom margins to ensure that the title and x axis title is visible.
                plt.subplots_adjust(bottom=0.15)
                plt.subplots_adjust(top=0.85)
                # Place the legend outside the axis
                plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                # Adjust figure size to increase width according to the outside legend width
                fig_size = plt.gcf().get_size_inches()
                fig_size[0] = fig_size[0] + fig_size[0]*(legend_width_scaled) # Increase width by the padding
                plt.gcf().set_size_inches(fig_size)
            else:
                plt.legend()
        # Use normalised species and orders lists to create the x-axis title with integral symbols.
        x_label="$\u03A3"
        for i in range(len(orders)):
            x_label=x_label+f"[{normalised_species[i]}]^{{{orders[i]}}}"
        x_label=x_label+"\u0394t$."
        # Define the x label. The size of the xlabel is scaled with "xaxis_size_scaler" if defined. 
        if xaxis_size_scaler!=None:
            plt.xlabel(x_label,fontsize=10*xaxis_size_scaler)
        else:
            plt.xlabel(x_label)
        # Define the y axis label using the output species name and the unit of its concentrations.
        y_label=f'[{output_species}] / {y_unit}'
        # Add shifted to the y axis label if all output species profiles don't start with the same value.
        if self.shifted_output==True:
            y_label='Shifted '+y_label
        plt.ylabel(y_label)
        # Define the plot title.
        if title==None:
            title_string=f'VTNA overlay plot for {output_species}'
        else:
            title_string=title
        plt.title(title_string)
        # Add gridlines if specified.
        if grid:
            plt.grid()
        # Set remaining plot parameters. 
        plt.xlim(0)
        plt.xticks(rotation=xtick_rotation) 
        plt.ylim(0)
        # Apply scientific x axis notation if selected.
        if xaxis_notation in ['Scientific','scientific','sci','Sci']:
            plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        elif xaxis_notation in ['Normal','normal']:
            plt.ticklabel_format(style='plain', axis='x')
        # Set tight layout if the legend is not placed outise the axis. 
        if not legend_outside:
            plt.tight_layout()
        if plt_show:
            plt.show()
        else:
            return fig
    
    def plot_VTNA_with_overlay_score(self,y_unit='M',size_scaler=1,title=None,grid=False,\
                                     xaxis_size_scaler=None,DP_scaler=1,fit_metric='RMSE',\
                                     deg=5,constraint='monotonic',show_fit_function=True,\
                                     xtick_rotation=0,show_overlay_score=True, show_legend=True,
                                     legend_outside=False,extra_legend=True,plt_show=True,
                                     linewidth=1,scale_xaxis=False,xaxis_notation='Automatic'):
        """
        Plots the VTNA overlay plot for the kinetic data according to the specifications of the \
        VTNA selection dictionary. Also plots the total fit used to generate the overlay score \
        using the specified fit metric, degree and constraint. 
        Args: 
            y_unit (str): String defining the output concentration unit for the overlay plots. 
            size_scaler (int or float): Factor by which the overlay plot size is adjusted. \
            title (str): String defining the title of the plot. 
            grid (bol): Set to True to create a grid based on x and y axis tick locations. 
            xaxis_size_scaler (int or float): Factor by which the size of the xaxis label is adjusted. 
            shifted_output (bol): Set to True to add "shifted" to the y axis title if a reactant has 
            been defined as the output reaction species. \
            fit_metric (str): Goodness-of-fit metric chosen to define the overlay score value. Either 'R2', 'RMSE' or 'SE'. 
            deg (int): Polynomial degree used in the global fitting prodecure to calculate the overlay score. 
            constraint (str): Defines the constraint of the global fit used by explore_orders() to calculate the overlay score. 
            show_fit_function (bol): Set to True to visualise the fit function underlying the overlay score. 
            show_overlay_score (bol): Set to True to include the overlay score under the title of the plot. 
            xtick_rotation (int or float): Determines the degree of rotation of the xaxis tick numbers. 
            show_legend (bol): Set to True to show the legend for different reaction profiles or to False \
            to omit it. 
            legend_outside (bol): Set to True to show the legend for the different reaction profiles outside \
            the main plot. Otherwise, it will be placed automatically at the most open space within the plot \
            automatically. 
            extra_legend (bol): Set to True to show the legend specifying the nature of the fit function. \
            Can be set to False to not show this extra legend. 
            plt_show (bol): Set to True to apply plt.show() at the end of the code to ensure that the plot is \
            generated. Otherwise if set to False, the figure can be modified externally after the execution of the code, giving \
            greater custumisability. 
            linewidth (float): Width of the line that connects the datapoints of concentration profiles. 
            scale_xaxis (bol): Set to True to scale the x axis values (normalised time axis) from 0 to 1. 
            xaxis_notation (string): Determines the notation of the x axis. Can be set to 'Normal' to avoid scientific \
            notation, 'Scientific' to ensure scienfific notation eg. 1e2 instead of 100, or 'Automatic' to get the default setting.
        """
        plt.rcParams['font.family'] = 'DejaVu Sans'
        # If neither the fit function nor the overlay score is set to be shown, use plot_VTNA instead.
        if not show_fit_function and not show_overlay_score:
            self.plot_VTNA(y_unit=y_unit,size_scaler=size_scaler,title=title,grid=grid,xaxis_size_scaler=xaxis_size_scaler,
            DP_scaler=DP_scaler,show_legend=show_legend,legend_outside=legend_outside,xtick_rotation=xtick_rotation,linewidth=linewidth,
            xaxis_notation=xaxis_notation)
            return
        # Determine the size of the overlay plot.
        fig,ax = plt.subplots(figsize=(5*size_scaler,4*size_scaler))
        # Determine the normalised and output reaction species from the selection dictionary.
        output_species=self.VTNA_selection['output_species']
        normalised_species=[]
        # Create a list of the order values of the different normalised reaction speices. 
        orders=[]
        for ns in self.VTNA_selection['normalised_species'].keys():
            normalised_species.append(ns)
            orders.append(str(self.VTNA_selection['normalised_species'][ns]))
        tables_scaled=data_scaling_y(self.result,self.VTNA_selection)
        # If scale_xaxis is set to true, scale the transformed time values by the largest across the dataset. 
        if scale_xaxis:
            tTs=[]
            for i in tables_scaled.keys():
                tTs=tTs+list(tables_scaled[i]['tT'])
            scaler=max(tTs)
            for i in tables_scaled.keys():
                tables_scaled[i]['tT']=np.array(tables_scaled[i]['tT'])/scaler
        # Use normalised species and orders lists to create the x-axis title with integral symbols.
        x_label="$\u03A3"
        for i in range(len(orders)):
            x_label=x_label+f"[{normalised_species[i]}]^{{{orders[i]}}}"
        x_label=x_label+"\u0394t$."
        # Define the x label. The size of the xlabel is scaled with "xaxis_size_scaler" if defined. 
        if xaxis_size_scaler!=None:
            ax.set_xlabel(x_label,fontsize=10*xaxis_size_scaler)
        else:
            ax.set_xlabel(x_label)
        # Set the datapoint size scaler to 1 if an incorrect value has been assigned.
        if DP_scaler==None or isinstance(DP_scaler, str):
            DP_scaler=1
        # Check difference between first and last output values to judge whether fit polynomial should 
        # be monotonically increasing or decreasing.
        delta_output=0
        for RM in tables_scaled.keys():
            if (tables_scaled[RM]['output_scaled'][0]-tables_scaled[RM]['output_scaled'].to_numpy()[-1])>0:
                delta_output=delta_output-1
            else: 
                delta_output=delta_output+1
        fit_direction='dec'
        if delta_output>0:
            fit_direction='inc'
        # Plot the time-normalised reaction profiles for each experiment dataset in tables_scaled
        Data_points_y=[]
        Data_points_t=[]
        # Find the largest tT value in the tables_scaled dataframe if scale_xaxis is set to True.
        for i in tables_scaled.keys():
            x_values=[float(j) for j in tables_scaled[i]['tT'].to_numpy()]
            Data_points_t=Data_points_t+x_values
            ys=[float(j) for j in tables_scaled[i]['output_scaled'].to_numpy()]
            Data_points_y=Data_points_y+ys
            ax.plot(x_values,ys,marker='o',label=i,markersize=6*DP_scaler,linewidth=linewidth)
        # Generate a legend if specified.
        if show_legend:
            if legend_outside:
                # Place the legend outside the axis
                legend_1=ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                ax.add_artist(legend_1)
                # Obtain the legend width in display units.
                legend_box = legend_1.get_window_extent()
                legend_width = legend_box.width
                # Obtain the figure width in display units.
                fig_box = ax.figure.get_window_extent()
                fig_width_display = fig_box.width
                # See by what fraction the figure width needs to be extended for outside legend.
                legend_width_scaled = legend_width / fig_width_display
                # Adjust figure size to increase width according to the outside legend width
                fig_size = ax.figure.get_size_inches()
                fig_size[0] = fig_size[0] + fig_size[0] * (legend_width_scaled)  # Increase width by the padding
                ax.figure.set_size_inches(fig_size)
                # Obtain the new figure width.
                fig_width_display = fig_box.width
                # Obtain the relative length of the legend compared to the figure.
                legend_width_scaled_2 = legend_width / fig_width_display
                # Scale the right argument according to the legend and new figure width. 
                fig.subplots_adjust(right=(1 - legend_width_scaled_2)*0.97)
                fig.subplots_adjust(bottom=0.15)
                fig.subplots_adjust(top=0.85)
                # Adjust top and bottom margins to ensure that the title and x axis title is visible.
            else:
                # Define a legend for the plot.
                legend_1=ax.legend()
                ax.add_artist(legend_1)
        # Define the y_label:
        y_label=f'Scaled [{output_species}]'
        # Add shifted to the y axis label if all output species profiles don't start with the same value.
        if self.shifted_output==True:
            y_label='Shifted '+y_label
        ax.set_ylabel(y_label)
        # Define a grid of transformed time values for the fit function.
        tT_2=np.linspace(0,1,100).reshape((-1,1))
        tT_2_not_scaled=np.linspace(0,max(Data_points_t),100).reshape((-1,1))
        if constraint=='monotonic':
            # Initiate the monotonic polynomial fitting method.
            polyestimator = PolynomRegressor(deg=deg, regularization = None, lam = 0)
            monotone_constraint = Constraints(monotonicity=fit_direction)
            # Perform the polynomial fit of "deg" order and evaluate fitted function at the relevant "tT" values       
            Data_points_t_2=np.array([i/max(Data_points_t) for i in Data_points_t]).reshape((-1,1))
            polyestimator.fit(Data_points_t_2, Data_points_y, loss = 'l2',verbose=False, constraints={0: monotone_constraint})    
            # Predict the actual datapoints.
            pred_mon = polyestimator.predict(Data_points_t_2)
            # Generate the fit y-values to plot the line. 
            pred_mon_fit = polyestimator.predict(tT_2)
            # Add the fit function plot to the axis if specified. 
            if show_fit_function:
                fit_plot,=ax.plot(tT_2*max(Data_points_t),pred_mon_fit,marker='',label=f'Monotonic {deg}.\norder polynomial fit',linestyle='dotted',linewidth=2,color='black')
            # Call the score_fit function to evalue the goodness-of-fit using the fit measuremetn "fit_metric"
            FQ=score_fit(Data_points_y,pred_mon,fit_metric)
        # If constraint is set to 'origin', perform total fit with a linear curve through origin.
        elif constraint=='origin':
            params = curve_fit(origin_line, Data_points_t, Data_points_y)
            a=params[0]
            # Generate predicted y values to evaluate the goodness of fit. 
            y_pred=a*Data_points_t
            FQ=score_fit(Data_points_y,y_pred,fit_metric)
            if show_fit_function:
                fit_plot,=ax.plot(tT_2_not_scaled,tT_2_not_scaled*a,marker='',label='linear fit \nvia origin',linestyle='dotted',linewidth=2,color='black')
        else:
            # Warn that ordinary polynomial fitting is conducted if an unknown argument is given. 
            if constraint not in [None,'None']:
                warnings.warn("Ordinary polynomial fitting applied due to unknown constraint input", UserWarning)
            # Perform the fit. 
            polyfit=np.polyfit(Data_points_t,Data_points_y,deg=deg)
            # Predict the y-values of every datapoint to determine the overlay score. 
            pred_y=np.polyval(polyfit,Data_points_t)
            FQ=score_fit(Data_points_y,pred_y,fit_metric)
            # Generate y values to plot the fit function. 
            pred_y_fit=np.polyval(polyfit,tT_2_not_scaled)
            if show_fit_function:
                fit_plot,=ax.plot(tT_2_not_scaled,pred_y_fit,marker='',label=f'{deg}. order\npolynomial fit',linestyle='dotted',linewidth=2,color='black')
        FQ_round=round(FQ,5)
        # Define the plot title.
        if title==None:
            title_string=f'VTNA overlay plot for {output_species}'
        else:
            title_string=title
        # Add the overlay score under the plot title if show_overlay_score is set to True. 
        GOF_string,slope_string='',''
        if show_overlay_score:
            GOF_string=f'\nGOF ({fit_metric}): {FQ_round}'
        if deg==1:
            if not scale_xaxis:
                tTs=Data_points_t
            Data_points_y_no_scaling=[]
            for i in tables_scaled.keys():
                ys=[float(j) for j in tables_scaled[i][output_species].to_numpy()]
                Data_points_y_no_scaling.extend(ys)
            if constraint=='origin':
                params, covariance = curve_fit(origin_line, tTs, Data_points_y_no_scaling)
            else:
                params, covariance = curve_fit(simple_line, tTs, Data_points_y_no_scaling)
            slope=params[0]
            print(slope)
            std_err_slope = np.sqrt(np.diag(covariance))[0]
            print(covariance)
            print(std_err_slope)
            # Determine a sensible number of decimals to include for the slope value. 
            rounding=5
            order_of_magnitude_slope=round(np.log10(abs(slope)))
            rounding+=(-order_of_magnitude_slope)
            slope_string=f'\nSlope (non-scaled tT and y): {round(slope,rounding)} ± {round(std_err_slope,rounding)}'
        plt.title(title_string+GOF_string+slope_string)
        # Generate the legend for the fit function. 
        if show_fit_function and extra_legend:
            ax.legend(handles=[fit_plot],loc='lower center')
        # Add gridlines if specified.
        if grid:
            plt.grid()
        ax.tick_params(axis='x', rotation=xtick_rotation)
        # Apply scientific x axis notation if selected.
        if xaxis_notation in ['Scientific','scientific','sci','Sci']:
            ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        elif xaxis_notation in ['Normal','normal']:
            plt.ticklabel_format(style='plain', axis='x')
        # Set tight layout if the legend is not placed outise the axis. 
        if not legend_outside:
            plt.tight_layout()
        # Adjust top and bottom margins to ensure that the title and x axis title is visible. 
        plt.subplots_adjust(bottom=0.15)
        top_margin=0.85 if deg>1 or not show_overlay_score else 0.82
        plt.subplots_adjust(top=top_margin)
        plt.xlim(0)
        plt.ylim(0)
        if plt_show:
            plt.show()