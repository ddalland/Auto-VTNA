import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from num2words import num2words
import re
import mplcursors
import math

def check_kinetic_data(kinetic_data,bold_and_red=True,print_report=True):
    """
    Checks whether the kinetic data loaded has the right format for auto_vtna and creates a report \
    string with relevant information and warnings. The function checks that the kinetic_data \
    variable is a dictionary with dataframe elements, that the column titles are consistent across \
    experiment sheets, that there are more than 2 columns per experiment dataframe and that all \
    concentration and time values are numerical. To load the kinetic data correctly, first \
    create an excel file with the concentration profiles from each experiment in separate sheets.\
    Each sheet should start with a time column with a suitable column title in cell A1, followed by concentration \
    values over time for different reaction species, each with column titles referring to the compound name \
    in the first row of sheet (B1, C1, D1 etc.): kinetic_data=pd.read_excel('file_name.xlsx',sheet_name=None). \
    Alternatively, the data from each experiment can be loaded as separate csv files using \
    pd.read_csv('filename.csv') and then collected as a dictionary manually. This takes more time.\
    """
    # Define functions used to check if for numerical values in lists.
    def all_elements_not_convertible_to_float(my_list):
        return all(not isinstance(element, (int, float)) and not (isinstance(element, str) 
        and element.replace('.', '', 1).isdigit()) for element in my_list)
    def check_numerical(values):
        return [isinstance(val, (int, float)) for val in values]
    def count_convertible_to_float(my_list):
        return sum(isinstance(element, (int, float)) or (isinstance(element, str) and 
        element.replace('.', '', 1).isdigit()) for element in my_list)
    # Define function to remove escape sequences from the result string. 
    def remove_escape_sequences(input_string):
        # Define the pattern for escape sequences
        escape_sequence_pattern = re.compile(r'\x1b\[[0-9;]*m')
        # Use sub() to replace all matches with an empty string
        cleaned_string = escape_sequence_pattern.sub('', input_string)
        return cleaned_string
    error=None
    tick_symbol='\u2713'
    # Begin defining the report string and its first section. 
    report='Kinetic data check report.\n--------------------------\n\033[1mCheck 1: Correct \
loading of kinetic data from Excel.\033[0m\n'
    # Ensure that kinetic_data is a dictionary, end the check if it isn't. 
    has_empty_dataframes = any(df.empty for df in kinetic_data.values())
    if type(kinetic_data)!=dict or has_empty_dataframes:
        report=report+"\033[91mOBS:\033[0m The kinetic data is not loaded correctly. It should \
be a dictionary containing the different experiment datasets as dataframes\nMake sure that there \
are no empty sheets and that the excel file contains one experiment dataset per\nsheet with column \
titles for time and concentration of reaction species and load the data like this: \
pd.read_excel('file_name.xlsx',sheet_name=None)\n"
        report=report+f'\033[1mConclusion: the kinetic_data variable is not ready for use with the\
auto_vtna package.\033[0m'
        # Remove escape sequences from the result string if bol_and_red set to False. 
        if bold_and_red==False:
            report=remove_escape_sequences(report)
        if print_report:
            print(report)
        return report
    # Show that it has been confirmed that kinetic_data is a dictionary. 
    report=report+f'  {tick_symbol} The kinetic data is loaded as a dictionary.\n'
    # Check that each key:value pair in kinetic_data contains a pandas dataframe object.
    type_elements=[]
    for i in kinetic_data.keys():
        type_elements.append(str(type(kinetic_data[i])))
    if all(i=="<class 'pandas.core.frame.DataFrame'>" for i in type_elements):
        report=report+f'  {tick_symbol} Each element associated a key in the kinetic data dictionary \
is a pandas dataframe object.\n'
    # If each key:value pair in kinetic_data isn't a pandas dataframe object, end the check.
    else: 
        report=report+"\033[91mOBS:\033[0m Some of the elements in the kinetic data dictionary are \
not pandas dataframe objects. This will need to be fixed.\nMake sure the excel file contains one \
experiment dataset per\n sheet with column titles for time and concentration of reaction species \
and load the data like this: pd.read_excel('file_name.xlsx',sheet_name=None)\n"
        report=report+f'\033[1mConclusion: the kinetic_data variable is not ready for use with \
the auto_vtna package.\033[0m'
        # Remove escape sequences from the result string if bol_and_red set to False. 
        if bold_and_red==False:
            report=remove_escape_sequences(report)
        if print_report:
            print(report)
        return report
    # List the experiment labels from each dataframe in the kinetic_data dictionary.
    no_experiments=len(list(kinetic_data.keys()))
    # Create a list containing lists of column titles from each dataframe. 
    column_titles=[]
    for exp in kinetic_data.keys():
        column_titles.append(kinetic_data[exp].columns.tolist())
    # Create a list containing the number of columns for each dataframe. 
    no_columns=[]
    columns=[]
    for exp in kinetic_data.keys():
        no_columns.append(len(kinetic_data[exp].columns))
        columns.append(list(kinetic_data[exp].columns))
    # Check if there is only one dataframe in kinetic_data (not VTNA compatible).
    if len(no_columns)==1:
        error=True
        report=report+'\033[91mOBS:\033[0m The kinetic data only contains data from one \
experiment. This is not enough for VTNA.\n'
    # Check if the number of columns in each dataframe is the same. 
    if all(i==no_columns[0] for i in no_columns):
        report=report+f'  {tick_symbol} The dataframe element(s) in the kinetic data \
dictionary all contain {no_columns[0]} columns.\n'
        # Check if every dataframe contains only one column (not VTNA compatible).
        if 1 in no_columns:
            report=report+'\033[91mOBS:\033[0m The dataframe elements in the kinetic \
data dictionary contains only one column.\
            This is insufficient for VTNA.\n'
            report=report+f'\033[1mConclusion: the kinetic data is not ready for use \
with the auto_vtna package.\033[0m'
            # Remove escape sequences from the result string if bol_and_red set to False. 
            if bold_and_red==False:
                report=remove_escape_sequences(report)
            if print_report:
                print(report)
            return report
    # Check if the column titles are the same across every kinetic experiment.
    report=report+'\033[1mCheck 2: Column title consistency.\033[0m\n'
    # Check if the columns from each dataframe are the same length. 
    if all(i==no_columns[0] for i in no_columns):
        column_titles = np.array(column_titles)
        # Check if all the reaction species concentration columns are the same.
        are_equal = np.all(column_titles == column_titles[0, :], axis=0)
        no_columns=len(column_titles[0])
        if all(are_equal):
            report=report+f'  {tick_symbol} All experiment datasets have the same column titles.\n'
        else:
            # Update the error flag.
            error=True
            # Flag each column title that is not consistent across dataframes. 
            for i in range(no_columns):
                if are_equal[i]==False:
                    report=report+f'\033[91mOBS:\033[0m The {num2words(i+1, to="ordinal_num")} \
column title is different across experiments: {column_titles[:,i]} \n'
    # If the number of columns is inconsistent across dataframes, update the error flag.
    else:
        error=True
        # Add an OBS statement to the report string to informa that the number of columns is 
        # inconsistent and list the column titles of each dataframe. 
        report=report+ f'\033[91mOBS:\033[0m'
        report=report+f' the number of columns in the dataframe objects in \
kinetic_data are not consistent.\nThe column titles are:\n'
        for i in range(len(column_titles)):
            report=report+f"  - {list(kinetic_data.keys())[i]}: {column_titles[i]}\n"
    # Check if all kinetic data dataframes contain only numerical values.
    counter=0
    report=report+'\033[1mCheck 3: Numerical datapoints.\033[0m \n'
    for exp in kinetic_data.keys():
        for column in list(kinetic_data[exp].columns):
            # Create a list which is True for all the values in the column that are numerical values.
            num_column=[isinstance(val, (int, float)) and not math.isnan(val) for val in list(kinetic_data[exp][column])]
            # If no values in the column are non-numerical and therefore all values in num_column True,
            # Go to the next column. 
            if all(num_column):
                continue
            else:
                # If there are non-numerical values in the column, update the error flag and 
                # add an OBS statement to the report string. 
                error=True
                counter=counter+1
                report=report+ '\033[91mOBS:\033[0m'
                report=report+f' Row(s) number {[index+2 for index,value in enumerate(num_column) if value==False]}\
 in column "{column}" in dataset "{exp}" contains values which are not numerical. Needs to be fixed.\n'.replace('[','').replace(']','')
    # If the counter is still 0, this means all dataframes are fully numerical. 
    if counter==0:
        report=report+f'  {tick_symbol} Every experiment dataset contains only numerical\
values.\n'
    counter=0
    report=report+'\033[1mCheck 4: Positive datapoints.\033[0m \n'
    for exp in kinetic_data.keys():
        for column in list(kinetic_data[exp].columns):
            # Create a list which is True for all the values in the column that are not negative numerical values.
            any_negative_values = list(kinetic_data[exp][column].apply(lambda s: pd.to_numeric(s, errors='coerce') < 0))
            # If no values in the column are negative and therefore all values in any_negative_values True,
            # Go to the next column. 
            if not any(any_negative_values):
                continue
            else:
                # If there are non-numerical values in the dataframe, update the error flag and 
                # add an OBS statement to the report string. 
                error=True
                counter=counter+1
                report=report+ '\033[91mOBS:\033[0m'
                report=report+f' Row(s) number {[index+2 for index,value in enumerate(any_negative_values) if value==True]}\
 in column "{column}" in dataset "{exp}" contains negative values. Needs to be fixed.\n'.replace('[','').replace(']','')
    # If the counter is still 0, this means all dataframes are fully numerical. 
    if counter==0:
        report=report+f'  {tick_symbol} No datasets contain negative time or concentration values.\n'
    # Check if the column titles are the same across every kinetic experiment.
    
    # Proceed to check whether the time values in each sheet are monotonic. 
    report=report+'\033[1mCheck 5: Time column monotonicity.\033[0m\n'
    counter=0
    for exp in kinetic_data.keys():
        time_array=kinetic_data[exp][kinetic_data[exp].columns[0]].to_numpy()
        if not np.all(time_array[1:] > time_array[:-1]):
            report=report+ f'\033[91mOBS:\033[0m The time values in dataset {exp} are not strictly monotonically increasing. Needs to be fixed.\n'
            error=True
            counter=counter+1
    if counter==0:
        report=report+f'  {tick_symbol} The time column (first column) of the dataset in every sheet is monotonically increasing.\n'
    # Check if the column titles are the same across every kinetic experiment.
    report=report+f'\033[1mCheck 6: Time and concentration columns.\033[0m \n'
    check_for_floats=0
    # Create a flattened column title list. 
    columns_flat=[]
    for i in columns:
        columns_flat=columns_flat+i
    # Check if all column titles are strings. 
    if all(type(i)==str for i in columns_flat):
        report=report+f'  {tick_symbol} All dataframe elements in the kinetic data dictionary \
has column titles which are strings.\n' 
    else: 
        for exp_no in range(len(list(kinetic_data.keys()))):
            # Ensure that all column titles aren't strings that can be converted to numerical values.
            if all_elements_not_convertible_to_float(columns[exp_no]):
                continue
            else:
                # If some column titles are numerical, add an OBS statement to the report string.
                check_for_floats=1
                report=report+f"\033[91mOBS:\033[0m {count_convertible_to_float(columns[exp_no])}\
of the column titles of dataframe {list(kinetic_data.keys())[exp_no]} are numerical values.\n"
                # Update error as numerical column titles reflect incorrect kinetic_data format.
                error=True
    # If no column title is a numerical value, add this to the report string. 
    if check_for_floats==0:
        report=report+f'  {tick_symbol} The column titles of every experiment dataset are all \
strings.\n'
    # Define the time column titles from each experiment dataframe. 
    time_column_titles=[i[0] for i in columns]
    # Check if the time column title is consistent across dataframes. 
    if all(i==time_column_titles[0] for i in time_column_titles):
        time_column_titles=time_column_titles[0]
    # Make the user aware of the title of the columns that will be interpreted as time values. 
    report=report+f"\033[91mNB:\033[0m Make sure that the time columns are the first columns in \
each dataframe in kinetic_data. Currently, the first column titles are: {time_column_titles}.\n"
    # If no errors have been flagged, generate a positive conclusion (kinetic_data is ready)
    if error==None:
        report=report+f'\033[1mConclusion:\033[0m\033[0m The kinetic data has been loaded \
correctly if the time axis label is correct.'
    # Generate a conclusion warning if one or more errors in kinetic_data errors have been detected. 
    else:
        report=report+f"\033[91m\033[1mConclusion:\033[0m The kinetic data has not been loaded \
correctly.Make sure that the kinetic data Excel spreadsheet is organised correctly and try again \
using kinetic_data=pd.read_excel('file_name.xlsx',sheet_name=None)."
    # Remove escape sequences from the result string if bol_and_red set to False. 
    if bold_and_red==False:
        report=remove_escape_sequences(report)
    if print_report:
        print(report)
    return report

def make_VTNA_selection(data,normalised_species='NS',output_species='OS'):
    """
    Creates a default VTNA selection dictionary structure with keys for each \
    experiment name as well as the normalised_species and output_species. The \
    dictionary can then be copied by the user and information about the desired \
    VTNA filled in such as:
    - RM: Sub dictionary containing each experiment dataset (RM stands for reaction \
    mixture) with an associated list of omitted reaction species by list index. 
    - normalised_species: Sub dictionary containing information about which reaction \
    species' concentration profiles will be used to normalise the time axis. Each \
    normalised species has an associated reaction order value for time normalisation. 
    - output_species: The reaction species whose concentration profile will define \
    y-axis. 
    The keyword arguments normalised_species and output_species can be defined to \
    generate a VTNA selection dictionary with the relevant keys.
    """
    selection_dict={}
    selection_dict['RM']={}
    for i in data.keys():
        selection_dict['RM'][i]={}
        selection_dict['RM'][i]['omissions']=None
        if i==normalised_species:
            selection_dict['normalised_species']={}
            selection_dict['normalised_species'][i]=1
        if i==output_species:
            selection_dict['output_species']=output_species
    if output_species=='OS':
        selection_dict['output_species']={}
        selection_dict['output_species']='OS'
    if normalised_species=='NS':
        selection_dict['normalised_species']={}
        selection_dict['normalised_species']['NS']=1
    return selection_dict

def initial_concs(data):
    """
    Creates a Pandas DataFrame showing the initial concentrations of every reaction species across each \
    experiment. Requires that each experiment DataFrame in "data" have the same column titles.
    """
    # Obtain a list of the experiment names
    experiments=list(data.keys())
    # Obtain a list containing the reaction species with concentration profiles in the kinetic dataset.
    reaction_species=data[experiments[0]].columns[1:]
    # Loop through each experiment and obtain the intial concentration of each reaction species.
    initial_conc_dict={'Reaction species':reaction_species}
    for j in experiments:
        initial_conc_dict[f'{j} init. conc.']=data[j].iloc[0,1:].to_numpy()
    return pd.DataFrame(initial_conc_dict)

def plot_data(data, ylim=None, y_unit='M', t_unit=None, fig_size_scaler=1,plot_mode='together',
    legend_outside=False, DP_scaler=1,linewidth=1):
    """
    Inputs the kinetic data "data" that has been imported from an excel file using pd.read_excel() \
    as a dictionary of DataFrames containing the concentration profiles for each different excess \
    experiment. Plots the concentration profiles of each each experiment in separate plots unless \
    plot_mode is set to 'together' or 'scrollable'.
    Args:
        data (dict): Kinetic data as a dictionary of Pandas dataframes.
        ylim (float, optional): Y axis upper cutoff value so that Y∈[0, ylim]. 
        y_unit (str, optional): Concentration unit to be included in the Y axis label. 
        t_unit (str, optional): Time unit to be included in the X axis label. 
        fig_size_scaler (float or int, optional): Scaler to adjust the size of the graphs by the \
        selected factor.
        plot_mode (str): Either set to "scrollable" to generate a single graph which updates \
        to the next or previous graph by clicking right or left arrow keys, or 'together' to \
        plot graphs for each experiment together. If plot_mode is set to something else, \
        graphs will pop-up one after the other. 
        legend_outside (bol): Set to False by default to show legend inside the axis of the plot. \
        Can be set to True to show the legend outside the axis (top right). 
        DP_scaler (float, int): Factor by which the datapoints of the plots is scaled from a standard \
        value of 6. 
        linewidth (float, int): The linewidth parameter used to generate plots. Is by default set to 1. 
    """
    # Make the plot mode 'separate' if there is only one experiment in the dataset. 
    if len(list(data.keys()))==1:
        plot_mode='separate'
    # Set the datapoint size scaler to 1 if an incorrect value has been assigned.
    if DP_scaler==None or isinstance(DP_scaler, str):
        DP_scaler=1
    if plot_mode=='together':
        addition=0
        if legend_outside:
            addition=1
        # Check how many rows and columns there should be based on the number of dataframes in data.
        n,m = len(list(data.keys()))+addition,3
        if n>6:
            m=4
        num_cols = min(m, n)
        num_rows = (n + num_cols - 1) // num_cols
        # Add an exception if there are 9 plot objects:
        if len(list(data.keys()))+addition==9:
            num_cols,num_rows=3,3
        # Create a single figure and subplots with specified number of rows and columns.
        fig, axs = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=(7 * 0.45*fig_size_scaler * num_cols, 5 * 0.6*num_rows * fig_size_scaler),constrained_layout=True)
        for i, (rm_name, rm_data) in enumerate(data.items()):
            # Check if there's only one row
            if num_rows == 1:  
                ax = axs[i % num_cols]  # Use modulo to handle indexing in a 1D array
            else:
                row, col = divmod(i, num_cols)
                ax = axs[row, col]
            time = rm_data[rm_data.columns[0]]
            ax.set_title(rm_name)
            # Plot kinetic profiles for each reaction species for a given reaction mixture data set
            for j in rm_data.columns[1:]:
                ax.plot(time, rm_data[j], label=j, linestyle='-', marker='o',markersize=6*DP_scaler,linewidth=linewidth)
            ax.set_xlabel(rm_data.columns[0])
            # Add a unit to the time axis title is specified. 
            if t_unit is not None:
                ax.set_xlabel(f'{rm_data.columns[0]} / {t_unit}')
            ax.set_ylabel(f'Concentration / {y_unit}')
            # Define y-axis range up to ylim if defined.
            if ylim is not None:
                ax.set_ylim(0, ylim)
            # Add a legend if legend_outise is not set to true:
            if not legend_outside:
                ax.legend()
        # Add an extra axis object if legend_outside is set to True. 
        if legend_outside:
            i=i+1
            if num_rows == 1:  # Check if there's only one row
                ax = axs[i % num_cols]  # Use modulo to handle indexing in a 1D array
                handles, labels = axs[0].get_legend_handles_labels()
            else:
                row, col = divmod(i, num_cols)
                ax = axs[row, col]
                handles, labels = axs[0,0].get_legend_handles_labels()
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
            ax.legend(handles, labels, loc='center')
        # Remove empty subplots if there are any
        for i in range(len(data)+addition, num_cols * num_rows):
            fig.delaxes(axs.flatten()[i])
        plt.show()
    elif plot_mode=='scrollable':
        # Define a class for creating scrollable plots
        class ScrollablePlot:
            def __init__(self, data_dict, ylim=None, y_unit='M', t_unit=None, 
                         fig_size_scaler=1,legend_outside=legend_outside,DP_scaler=1,linewidth=1):
                # Extract data and labels from the given dictionary
                data_1 = [data_dict[i] for i, j in data_dict.items()]
                data_labels = list(data_dict.keys())
                # Initialize instance variables
                self.data = data_1
                self.labels = data_labels
                self.current_plot_index = 0
                self.ylim = ylim
                self.y_unit = y_unit
                self.t_unit = t_unit
                self.fig_size_scaler = fig_size_scaler
                self.legend_outside=legend_outside
                self.DP_scaler=DP_scaler
                self.linewidth=linewidth
                # Create a Matplotlib figure and axis
                self.fig, self.ax = plt.subplots()
                # Create a note to inform user of scrolling feature.
                note_text = 'NB: Use left or right arrows to scroll to other experiments.'
                self.fig.text(0.5, 0.015, note_text, ha='center', color='gray')
                # Plot the initial data
                self.plot_current_data()
                # Enable hovering using mplcursors
                mplcursors.cursor(hover=True)
                # Connect the key press event to the corresponding method
                self.fig.canvas.mpl_connect('key_press_event', self.on_key_press)

            def plot_current_data(self):
                # Clear the current axis
                self.ax.clear()
                # Get the data and label for the current plot
                current_data = self.data[self.current_plot_index]
                current_label = self.labels[self.current_plot_index]
                # Extract time and concentration data
                t = current_data.columns[0]
                time_column = current_data[t]
                # Plot each concentration profile
                for i in current_data.keys():
                    if i == t:
                        continue
                    self.ax.plot(time_column, current_data[i], label=i, marker='o',
                                 markersize=6*self.DP_scaler,linewidth=self.linewidth)
                # Set the figure size
                self.fig.set_size_inches(7 *0.8* self.fig_size_scaler, 5 *0.8* self.fig_size_scaler)
                # Set plot title
                self.ax.set_title(f'Concentration profiles in {current_label}')
                # Set x-axis label
                if self.t_unit is not None:
                    self.ax.set_xlabel(f'Time / {self.t_unit}')
                else:
                    self.ax.set_xlabel(t)
                # Set y-axis label
                self.ax.set_ylabel(f'Concentration / {self.y_unit}')
                # Set y-axis limit if specified
                if self.ylim is not None and isinstance(self.ylim, (float, int)):
                    self.ax.set_ylim(0, self.ylim)
                # Display legend
                if self.legend_outside:
                    # Calculate the width of the legend dynamically based on the longest text element
                    legend = self.ax.legend()
                    # Obtain the legend width in display units. 
                    legend_bbox = legend.get_window_extent()
                    legend_width = legend_bbox.width
                    # Obtain the figure width in display units.
                    fig_box = plt.gcf().get_window_extent()
                    fig_width_display = fig_box.width
                    # See by what fraction the figure width needs to be extended for outside legend.
                    legend_width_scaled=legend_width/fig_width_display
                    # Scale the right argument according to the legend width
                    self.fig.subplots_adjust(right=1 - legend_width_scaled)
                    # Place the legend outside the axis
                    self.ax.legend(loc='center left',bbox_to_anchor=(1, 0.8))
                    # Adjust figure size to increase width according to the outside legend width
                    fig_size = plt.gcf().get_size_inches()
                    fig_size[0] = fig_size[0] + fig_size[0]*(legend_width_scaled) # Increase width by the padding
                    plt.gcf().set_size_inches(fig_size)
                else:
                    self.ax.legend()
                # Adjust layout to prevent overlapping of titles and labels
                plt.tight_layout()
                # Redraw the plot
                plt.draw()
            def on_key_press(self, event):
                # Handle key presses for navigating between plots
                if event.key == 'right':
                    self.next_plot()
                elif event.key == 'left':
                    self.previous_plot()
            def next_plot(self):
                # Move to the next plot and update the display
                self.current_plot_index = (self.current_plot_index + 1) % len(self.data)
                self.plot_current_data()
            def previous_plot(self):
                # Move to the previous plot and update the display
                self.current_plot_index = (self.current_plot_index - 1) % len(self.data)
                self.plot_current_data()
        scrollable_plot = ScrollablePlot(data, ylim=ylim, y_unit=y_unit, t_unit=t_unit, 
        fig_size_scaler=fig_size_scaler,legend_outside=legend_outside,DP_scaler=DP_scaler,linewidth=linewidth)
        plt.show()
    else:
            # Identify the first key in the data dictionary
        RM1=list(data.keys())[0]
        # Identify the time column unit for the first reaction mixture data set
        t=data[RM1].columns[0]
        for i in data.keys():
            plt.figure(figsize=(7*0.8*fig_size_scaler, 5*0.8*fig_size_scaler))
            time=data[i][t]
            # plot kinetic profiles for each reaction species for a given reaction mixture data set
            for j in data[i].keys():
                if j==t:
                    continue
                plt.plot(time,data[i][j],label=j,linestyle='-',marker='o',markersize=6*DP_scaler,linewidth=linewidth)
                plt.title(i)
                plt.xlabel(t)
                if t_unit!=None:
                    plt.xlabel(f'{t} / {t_unit}')
                plt.ylabel(f'Concentration / {y_unit}')
            # Define y-axis range up to ylim if defined. 
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
                # Place the legend outside the axis
                plt.legend(loc='center left', bbox_to_anchor=(1, 0.8))
                # Adjust figure size to increase width according to the outside legend width
                fig_size = plt.gcf().get_size_inches()
                fig_size[0] = fig_size[0] + fig_size[0]*(legend_width_scaled) # Increase width by the padding
                plt.gcf().set_size_inches(fig_size)
            else:
                plt.legend()
            if ylim!=None:
                plt.ylim(0,ylim)
            #Use plt.show to complete plot and move to the next reaction mixture data set
            plt.tight_layout()
            plt.show()

def plot_data_together(data,species,ylim=None,y_unit='mM',fig_size_scaler=1,t_unit=None,
    legend_outside=False,DP_scaler=1,linewidth=1):
    """
    Inputs the kinetic data "data" that has been imported from an excel file using Pd.read_excel() \
    as a dictionary of DataFrames containing the concentration profiles for each different excess \
    experiment. Creates a plot with the concentration profiles of the selected reaction species from \
    each experiment. 
    Args:
        data (dict): Kinetic data as a dictionary of Pandas dataframes.
        species (list): List containing the reaction species whose concentration profiles \
        will be plotted. 
        ylim (float or int, optional): Y axis upper cutoff value so that Y∈[0,ylim]. 
        y_unit (str, optional): Concentration unit to be included in the Y axis label.
        fig_size_scaler (float or int, optional): Scaler to adjust the size of the graphs by the \
        selected factor.
        t_unit (str): Unit of the time axis. Gives a "/ " followed by the string as an addition to \
        the column title of the time column in the inputted kinetic dataset. 
        legend_outside (bol): Set to False by default to show legend inside the axis of the plot. \
        Can be set to True to show the legend outside the axis (top right). 
        DP_scaler (float, int): Factor by which the datapoints of the plots is scaled from a standard \
        value of 6. 
        linewidth (float, int): The linewidth parameter used to generate plots. Is by default set to 1. 
    """
    # Identify the first key in the data dictionary
    RM1=list(data.keys())[0]
    # Identify the time column unit for the first reaction mixture data set
    t=data[RM1].columns[0]
    # Set the datapoint size scaler to 1 if an incorrect value has been assigned.
    if DP_scaler==None or isinstance(DP_scaler, str):
        DP_scaler=1
    t_label=t
    if t_unit!=None:
        t_label=t_label+f' / {t_unit}'
    fig, ax = plt.subplots(figsize=(7*0.9*fig_size_scaler, 5*0.9*fig_size_scaler))
    for i in data.keys():
        time=data[i][t]
        # plot kinetic profiles for the chosen species for each reaction mixture data set
        for j in data[i].keys():
            if j==t or j not in species:
                continue
            ax.plot(time,data[i][j],label=f'{i}, {species}',linestyle='-',marker='o',
                    markersize=6*DP_scaler,linewidth=linewidth)
            ax.set_title(f'concentration profiles for {species}')
            ax.set_xlabel(t_label)
            ax.set_ylabel(f'concentration / {y_unit}')
        # Define y-axis range up to ylim if defined.
        if ylim!=None:
            ax.set_ylim(0,ylim)
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
        # Place the legend outside the axis
        plt.legend(loc='center left',bbox_to_anchor=(1, 0.8))
        # Adjust figure size to increase width according to the outside legend width
        fig_size = plt.gcf().get_size_inches()
        fig_size[0] = fig_size[0] + fig_size[0]*(legend_width_scaled) # Increase width by the padding
        plt.gcf().set_size_inches(fig_size)
    else:
        ax.legend()
    plt.tight_layout()
    plt.show()

def plot_data_MB(data, species, stochiometry, scaler=1, ylim=None, t_unit=None, y_unit='mM', fig_size_scaler=1,
                 plot_mode='together',legend_outside=False,DP_scaler=1,linewidth=1):
    """
    Inputs the kinetic data "data" that has been imported from an excel file using Pd.read_excel() \
    as a dictionary of DataFrames containing the concentration profiles for each different excess \
    experiment. Plots the concentration profiles of each each experiment with the mass balance relative \
    to time=0 of the chosen reaction species and their stochiometries. The graphs can be plotted in one \
    figure (plot_mode = 'together'), in separate pop-up figures (plot_mode='separate') or in a single \
    figure that can be updated upon 
    Args:
        data (dict): Kinetic data as a dictionary of Pandas dataframes.
        species (list): List containing the reaction species whose concentration profiles \
        will be used to calculate the mass balance. 
        stociometry (list): List containing the stociometries corresponding to the species \
        list.
        ylim (float or int, optional): Y axis upper cutoff value so that Y∈[0,ylim]. 
        y_unit (str, optional): Concentration unit to be included in the Y axis label.
        fig_size_scaler (float or int, optional): Scaler to adjust the size of the graphs by the \
        selected factor.
        plot_mode (str): Either set to "scrollable" to generate a single graph which updates \
        to the next or previous graph by clicking right or left arrow keys, or 'together' to \
        plot graphs for each experiment together. If plot_mode is set to something else, \
        graphs will pop-up one after the other. 
        legend_outside (bol): Set to False by default to show legend inside the axis of the plot. \
        Can be set to True to show the legend outside the axis (top right). 
        DP_scaler (float, int): Factor by which the datapoints of the plots is scaled from a standard \
        value of 6. 
        linewidth (float, int): The linewidth parameter used to generate plots. Is by default set to 1. 
    """
    # Make the plot mode 'separate' if there is only one experiment in the dataset. 
    if len(list(data.keys()))==1:
        plot_mode='separate'
    # Find the title of the time column of the experiment datasets
    RM1 = list(data.keys())[0]
    t = data[RM1].columns[0]
    # Set the datapoint size scaler to 1 if an incorrect value has been assigned.
    if DP_scaler==None or isinstance(DP_scaler, str):
        DP_scaler=1
    MB_title_string = ''
    # Define the title of the mass balance axis
    for mol in range(len(species)):
        if mol != 0 and species[mol] == species[-1]:
            MB_title_string = MB_title_string + 'and '
        MB_title_string = MB_title_string + f'{species[mol]}, '
    if len(species) == 1:
        MB_title_string = f'{species[0]}. '
    if plot_mode=='together':
        addition=0
        if legend_outside:
            addition=1
        # Check how many rows and columns there should be based on the number of dataframes in data.
        n,m = len(list(data.keys()))+addition,3
        if n>6:
            m=4
        num_cols = min(m, n)
        num_rows = (n + num_cols - 1) // num_cols
        # Add an exception if there are 9 plot objects:
        if len(list(data.keys()))+addition==9:
            num_cols,num_rows=3,3
        # Create a single figure and subplots with specified number of rows and columns.
        fig, axes = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=(7.5 * 0.45*fig_size_scaler * num_cols, 4.5 * 0.6*num_rows * fig_size_scaler),constrained_layout=True)
        axes = axes.flatten() 
        for idx, (i, ax) in enumerate(zip(data.keys(), axes)):
            Concs = []
            largest_conc = []
            time = data[i][t]
            # Create a curve for each reaction species
            for j in data[i].keys():
                if j == t:
                    continue
                ax.plot(time, data[i][j], label=j, linestyle='-', marker='o', markersize=6*DP_scaler,linewidth=linewidth)
                ax.set_title(i)
                ax.set_xlabel(t)
                ax.set_ylabel(f'concentration / {y_unit}')
                largest_conc.append(max(data[i][j].to_numpy()))
                if t_unit is not None:
                    ax.set_xlabel(f'{t} / {t_unit}')
                if j == t or j not in species:
                    continue
                # Add concentration*stochiometry array for calculation of mass balance over time
                Concs.append(data[i][j].to_numpy() * stochiometry[species.index(j)])
            # Use Concs to create mass balance sum for each time point
            sums = np.zeros(len(Concs[0]))
            for i in Concs:
                sums = sums + i
            # Use sums to find the % mass balance at each time relative to the initial value.
            # If the first sum is 0, use the final sum instead. (could occur if only product(s) are chosen)
            if sums[0]>0:
                MB = sums / sums[0] * 100 * scaler
            else:
                MB = sums / sums[-1] * 100 * scaler
            ax2 = ax.twinx()
            ax2.set_ylabel(f'Mass balance for:\n {MB_title_string[:-2]}', color="purple")
            ax2.plot(time, MB, label="Mass balance %", linestyle='-', marker='o', color="purple",markersize=6*DP_scaler)
            ax2.plot([0, max(time)], [100, 100], linestyle='--', marker='', color='purple')
            ax2.legend(loc='upper center')
            # Ensure appropriate y limits with space for legends.
            ax2.set_ylim(0, max(MB) + 28)
            if ylim is not None and type(ylim) in [float,int]:
                ax.set_ylim(0, ylim)
            else:
                ax.set_ylim(0, max(largest_conc) + max(largest_conc) * 0.15)
            if not legend_outside:
                ax.legend()
        ax = axes[idx+1]
        handles, labels = axes[0].get_legend_handles_labels()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.legend(handles, labels, loc='center left')
        for i in range(len(data)+addition, num_cols * num_rows):
          fig.delaxes(axes.flatten()[i])
        # Adjust layout to prevent overlap
        plt.tight_layout()
        plt.show()
    elif plot_mode=='scrollable':
        class ScrollablePlot:
            def __init__(self, data_dict,species, stochiometry,y_unit='M',t_unit=None,fig_size_scaler=1,
                         ylim=None,scaler=1,legend_outside=False,DP_scaler=1,linewidth=1):
                data_1=[]
                for i,j in data_dict.items():
                    data_1.append(data_dict[i])
                data_labels=list(data_dict.keys())
                self.data = data_1
                self.labels = data_labels
                self.current_plot_index = 0
                self.y_unit=y_unit
                self.t_unit=t_unit
                self.ylim=ylim
                self.fig_size_scaler=fig_size_scaler
                self.species=species
                self.stochiometry=stochiometry
                self.MB_title_string=MB_title_string
                self.legend_outside=legend_outside
                self.DP_scaler=DP_scaler
                self.linewidth=linewidth
                mplcursors.cursor(hover=True)
                self.fig, self.ax = plt.subplots(figsize=(7 *0.8* fig_size_scaler, 5 *0.8* fig_size_scaler))
                # Create a note to inform user of scrolling feature.
                note_text = 'NB: Use left or right arrows to scroll to other experiments.'
                self.fig.text(0.5, 0.015, note_text, ha='center', color='gray')
                self.ax2=self.ax.twinx()
                self.plot_current_data()
                self.fig.canvas.mpl_connect('key_press_event', self.on_key_press)

            def plot_current_data(self):
                self.ax.clear()
                self.ax2.clear()
                current_data = self.data[self.current_plot_index]
                current_label = self.labels[self.current_plot_index]
                t=current_data.columns[0]
                time_column=current_data[t]
                Concs=[]
                largest_conc=[]
                for j in current_data.keys():
                    if j==t:
                        continue
                    self.ax.plot(time_column,current_data[j],label=j,linestyle='-',marker='o',
                                 markersize=6*self.DP_scaler,linewidth=linewidth)
                    self.ax.set_title(f'Concentration profiles in {current_label}')
                    self.fig.set_size_inches(7*0.8*self.fig_size_scaler, 5*0.8*self.fig_size_scaler)
                    self.ax.set_xlabel(t)
                    self.ax.set_ylabel(f'concentration / {self.y_unit}')
                    largest_conc.append(max(current_data[j].to_numpy()))
                    if self.t_unit!=None:
                        self.ax.set_xlabel(f'Time / {self.t_unit}')
                    if j==t or j not in self.species:
                        continue
                    # Add concentration*stochiometry array for calculation of mass balance over time
                    Concs.append(current_data[j].to_numpy()*self.stochiometry[self.species.index(j)])
                # Use Concs to create mass balance sum for each time point
                sums=np.zeros(len(Concs[0]))
                for i in Concs:
                    sums=sums+i
                # Use sums to find the % mass balance at each time relative to the intial value.
                MB=sums/sums[0]*100* scaler
                self.ax2.set_ylabel(f'Mass balance for:\n {self.MB_title_string[:-2]}',color="purple")
                self.ax2.yaxis.set_label_position('right')
                self.ax2.plot(time_column,MB,label="Mass balance %",linestyle='-',marker='o',color="purple",markersize=6*self.DP_scaler)
                self.ax2.plot([0,max(time_column)],[100,100],linestyle='--',marker='',color='purple')
                self.ax2.legend(loc='upper center')
                # Ensure appropriate y limits with space for legends. 
                self.ax2.set_ylim(0,max(MB)+28)
                if self.ylim!=None and type(self.ylim) in [float,int]:
                    self.ax.set_ylim(0,self.ylim)
                else:
                    self.ax.set_ylim(0,max(largest_conc)+max(largest_conc)*0.15)
                if self.legend_outside:
                    # Calculate the width of the legend dynamically based on the longest text element
                    legend = self.ax.legend()
                    # Obtain the legend width in display units. 
                    legend_bbox = legend.get_window_extent()
                    legend_width = legend_bbox.width
                    # Obtain the figure width in display units.
                    fig_box = plt.gcf().get_window_extent()
                    fig_width_display = fig_box.width
                    # See by what fraction the figure width needs to be extended for outside legend.
                    legend_width_scaled=legend_width/fig_width_display
                    # Scale the right argument according to the legend width
                    self.fig.subplots_adjust(right=1 - legend_width_scaled-0.07)
                    # Place the legend outside the axis
                    legend_position_scaler=fig_size_scaler**0.3 if fig_size_scaler<1 else 1
                    self.ax.legend(loc='center left',bbox_to_anchor=(1.14/legend_position_scaler, 0.8))
                    # Adjust figure size to increase width according to the outside legend width
                    fig_size = plt.gcf().get_size_inches()
                    fig_size[0] = fig_size[0] + fig_size[0]*(legend_width_scaled) # Increase width by the padding
                    plt.gcf().set_size_inches(fig_size)
                else:
                    self.ax.legend()
                plt.tight_layout()
                plt.draw()

            def on_key_press(self, event):
                if event.key == 'right':
                    self.next_plot()
                elif event.key == 'left':
                    self.previous_plot()

            def next_plot(self):
                self.current_plot_index = (self.current_plot_index + 1) % len(self.data)
                self.plot_current_data()

            def previous_plot(self):
                self.current_plot_index = (self.current_plot_index - 1) % len(self.data)
                self.plot_current_data()
        scrollable_plot = ScrollablePlot(data,species,stochiometry,fig_size_scaler=fig_size_scaler,
        t_unit=t_unit,y_unit=y_unit,scaler=scaler,ylim=ylim,legend_outside=legend_outside,linewidth=linewidth, DP_scaler=DP_scaler)
        plt.show()
    else: 
        # Create a twin plot for each experiment dataset
        for i in data.keys():
            Concs=[]
            largest_conc=[]
            fig,ax = plt.subplots(figsize=(7*0.8*fig_size_scaler, 5*0.8*fig_size_scaler))
            fig_size = plt.gcf().get_size_inches()
            time=data[i][t]
            # Create a curve for each reaction species
            for j in data[i].keys():
                if j==t:
                    continue
                ax.plot(time,data[i][j],label=j,linestyle='-',marker='o',
                        markersize=6*DP_scaler,linewidth=linewidth)
                ax.set_title(i)
                ax.set_xlabel(t)
                ax.set_ylabel(f'concentration / {y_unit}')
                largest_conc.append(max(data[i][j].to_numpy()))
                if t_unit!=None:
                    plt.xlabel(f'{t} / {t_unit}')
                if j==t or j not in species:
                    continue
                # Add concentration*stochiometry array for calculation of mass balance over time
                Concs.append(data[i][j].to_numpy()*stochiometry[species.index(j)])
            if legend_outside:
                # Calculate the width of the legend dynamically based on the longest text element
                legend = ax.legend()
                # Obtain the legend width in display units. 
                legend_bbox = legend.get_window_extent()
                legend_width = legend_bbox.width
                # Obtain the figure width in display units.
                fig_box = plt.gcf().get_window_extent()
                fig_width_display = fig_box.width
                # See by what fraction the figure width needs to be extended for outside legend.
                legend_width_scaled=legend_width/fig_width_display
                # Scale the right argument according to the legend width
                fig.subplots_adjust(right=1 - legend_width_scaled-0.07)
                # Place the legend outside the axis
                legend_position_scaler=fig_size_scaler**0.3 if fig_size_scaler<1 else 1
                ax.legend(loc='center left',bbox_to_anchor=(1.14/legend_position_scaler, 0.8))
                # Adjust figure size to increase width according to the outside legend width
                fig_size = plt.gcf().get_size_inches()
                fig_size[0] = fig_size[0] + fig_size[0]*(legend_width_scaled) # Increase width by the padding
                plt.gcf().set_size_inches(fig_size)
            else:
                ax.legend()
            # Use Concs to create mass balance sum for each time point
            sums=np.zeros(len(Concs[0]))
            for i in Concs:
                sums=sums+i
            # Use sums to find the % mass balance at each time relative to the intial value.
            MB=sums/sums[0]*100*scaler
            ax2=ax.twinx()
            ax2.set_ylabel(f'Mass balance for:\n {MB_title_string[:-2]}',color="purple")
            ax2.plot(time,MB,label="Mass balance %",linestyle='-',marker='o',color="purple",markersize=6*DP_scaler)
            ax2.plot([0,max(time)],[100,100],linestyle='--',marker='',color='purple')
            ax2.legend(loc='upper center')
            # Ensure appropriate y limits with space for legends. 
            ax2.set_ylim(0,max(MB)+28)
            if ylim!=None and type(ylim) in [float,int]:
                ax.set_ylim(0,ylim)
            else:
                ax.set_ylim(0,max(largest_conc)+max(largest_conc)*0.15)
            fig.subplots_adjust(right=0.85)
            plt.show()

