from auto_vtna.VTNA_functions import same_excess
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import copy

class Same_Excess:
    def __init__(self, data, selected_experiments, selected_species):
        # Initialize with data and perform Same_Excess calculations
        self.data = data
        self.selected_experiments = selected_experiments
        self.selected_species = selected_species
        # Calculate the shifted time axis.
        results=same_excess(
            self.data, self.selected_experiments, self.selected_species)
        if isinstance(results, str):
            self.result=None
            self.error=results
            return
        else:
            self.result, self.low_initial_conc_experiment = results
            self.error=None
            
    def plot_same_excess(self, y_unit='M', size_scaler=1, title=None, grid=False,
                         xaxis_size_scaler=None, DP_scaler=1, xtick_rotation=0,
                         show_legend=True, legend_outside=False, plt_show=True,
                         linewidth=1, xaxis_notation='Automatic',
                         title_fontsize=10):
        """
        Plots the Same excess plot for two selected experiments to see whether product inhibition or catalyst decomposition is present.
        """
        plt.rcParams['font.family'] = 'DejaVu Sans'
        # Check if the same excess calculation worked:
        if isinstance(self.result, str):
            print(self.result)
            return
        # Determine the shape of the overlay plot
        fig_size = (5 * size_scaler, 4 * size_scaler)
        fig, ax = plt.subplots(figsize=fig_size)
        
        # Set the datapoint size scaler to 1 if an incorrect value has been assigned.
        if DP_scaler is None or isinstance(DP_scaler, str):
            DP_scaler = 1
        
        # Define the time label.
        time_label = self.result[self.selected_experiments[0]].columns[0]

        # Identify initial shifted time value:
        initial_t_shift=self.result[self.low_initial_conc_experiment][time_label][0]
        # Create a copy of the data with the low initial concentration profile shifted back:
        data2 = copy.deepcopy(self.result)

        # Access the specific DataFrame in the dictionary
        df = data2[self.low_initial_conc_experiment]

        # Modify the time_label column by subtracting the initial_t_shift value
        df[time_label] = df[time_label] - initial_t_shift

        # Assign the modified DataFrame back to the dictionary (optional, since df is already a reference)
        data2[self.low_initial_conc_experiment] = df

        # Plot the same excess and standard experiments with a shifted time axis for the former.
        plots = {}
        for exp in self.result.keys():
            if exp == self.low_initial_conc_experiment:
                label = f"{exp} (Shifted time axis)"
            else:
                label = exp
            x_values = [float(j) for j in self.result[exp][time_label].to_numpy()]
            ys = [float(j) for j in self.result[exp][self.selected_species].to_numpy()]
            plot, = ax.plot(x_values, ys, marker='o', label=label, markersize=6 * DP_scaler, linewidth=linewidth)
            plots[exp] = plot
        
        # Set up the slider for horizontal shifting
        ax_shift = plt.axes([0.3, 0.01, 0.6, 0.03], facecolor='lightgoldenrodyellow')
        slider = Slider(ax=ax_shift, label='Shift S.E. profile', valmin=0, valmax=initial_t_shift*4, valinit=initial_t_shift, valstep=initial_t_shift/1000)
        
        def update(val):
            shift = slider.val
            x_values_shifted = [float(j) + shift for j in data2[self.low_initial_conc_experiment][time_label].to_numpy()]
            plots[self.low_initial_conc_experiment].set_xdata(x_values_shifted)
            fig.canvas.draw_idle()
        
        slider.on_changed(update)
        
        # Adjust layout to accommodate the slider
        plt.subplots_adjust(bottom=0.18)
        
        if show_legend:
            if legend_outside:
                # Place the legend outside the axis
                ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            else:
                ax.legend()
        
        # Define the x and y axis labels
        ax.set_xlabel(time_label)
        y_label = f'[{self.selected_species}] / {y_unit}'
        ax.set_ylabel(y_label)
        
        # Define the plot title
        title_string = title if title else f'Same excess plot for {self.selected_species}'
        ax.set_title(title_string, fontsize=title_fontsize)
        
        # Add gridlines if specified
        if grid:
            ax.grid()
        
        # Apply x-axis notation if specified
        if xaxis_notation.lower() in ['scientific', 'sci']:
            ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        elif xaxis_notation.lower() in ['normal']:
            ax.ticklabel_format(style='plain', axis='x')
        
        # Add horizontal and vertical lines at y=0 and x=0
        ax.axhline(0, color='black', linestyle='-',linewidth=0.7)
        ax.axvline(0, color='black', linestyle='-',linewidth=0.7)
        
        if plt_show:
            plt.show()
        else:
            return fig
