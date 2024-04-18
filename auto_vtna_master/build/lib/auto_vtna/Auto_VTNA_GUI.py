import PySimpleGUI as sg
import threading
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import auto_vtna
import tkinter as tk
from pandastable import Table
from auto_vtna.Normal_VTNA import Normal_VTNA
from auto_vtna.VTNA_functions import VTNA_omissions
from auto_vtna.Automatic_VTNA import Automatic_VTNA
from openpyxl.drawing.image import Image
import copy
import re
import csv
import tkinter as tk
from tkinter import scrolledtext

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

def overlay_plot(kinetic_data,sheets,NS_dict,output_reaction_species,tt_scaler=1,y_unit='M',\
                 size_scaler=1,grid=False,title=None,fit_metric='RMSE',constraint='monotonic',\
                 deg=5,show_fit_function=False,show_overlay_score=False,DP_scaler=1,xtick_rotation=0,\
                 show_legend=True,legend_outside=False):
    ''' Creates a VTNA_selection dictionary for the GUI settings stored in the arguments
    "sheets", "NS_dict" and "output_reaction_species", calculates the transformed time axis 
    for the kinetic data "kinetic_data" and generates the overlay plot for the selected plot settings.
    '''
    # Generate a VTNA selection dictionary template.
    Default_VTNA_selection_dictionary=auto_vtna.make_VTNA_selection(kinetic_data)
    VTNA_selection_dictionary={}
    # Generate a sub-dictionary to specify which experiments to include based on sheet selections.
    VTNA_selection_dictionary['RM']={}
    for i in sheets:
        VTNA_selection_dictionary['RM'][i]=Default_VTNA_selection_dictionary['RM'][i]
    # Update the normalised species and output species for generating the overlay plot based on GUI selections.
    VTNA_selection_dictionary['normalised_species']=NS_dict
    VTNA_selection_dictionary['output_species']=output_reaction_species
    # Check if the output species selected will require concentration profiles to be shifted vertically. 
    shifted_output=False
    if kinetic_data[sheets[0]][output_reaction_species].to_numpy()[0]>kinetic_data[sheets[0]][output_reaction_species].to_numpy()[-1]:
        shifted_output=True
    # Calculate the normalised time axis using Normal_VTNA from auto_vtna.
    Run_VTNA=Normal_VTNA(kinetic_data,VTNA_selection_dictionary)
    if constraint=='through origin':
        constraint='origin'
    # Generate the overlay plot.
    Run_VTNA.plot_VTNA_with_overlay_score(y_unit=y_unit,size_scaler=size_scaler,grid=grid, \
                  xaxis_size_scaler=tt_scaler,title=title,fit_metric=fit_metric,constraint=constraint,\
                  deg=deg,show_fit_function=show_fit_function,show_overlay_score=show_overlay_score,\
                    DP_scaler=DP_scaler,xtick_rotation=xtick_rotation,show_legend=show_legend,
                    legend_outside=legend_outside)
    
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

def run_calculation_in_thread(window, data, selected_sheet, selected_columns, selected_result_columns, calculation_settings):
    # NB: this function was created and annotated using OpenAI's ChatGTP. 
    # Create a thread for the calculation
    calculation_thread = CalculationThread(target=run_calculation, args=(data, selected_sheet, selected_columns, selected_result_columns, calculation_settings))
    calculation_thread.start()
    # Update the GUI to indicate that the calculation is underway
    popup_layout = [[sg.Text("Calculation is underway. Please wait...")]]
    popup_window = sg.Window("Please Wait", popup_layout, finalize=True)
    while True:
        # Read events from the popup window with a timeout of 100 milliseconds
        event, values = popup_window.read(timeout=100)  
        # Check if the popup window is closed or the calculation thread has finished
        if event == sg.WIN_CLOSED or not calculation_thread.is_alive():
            break
    # Close the popup window
    popup_window.close()
    # Wait for the calculation thread to complete
    calculation_thread.join()
    # Retrieve the return values from the calculation thread
    calculation_result = calculation_thread.result
    return calculation_result

def plot_kinetic_data(data, stochiometry_list, ylim=None, y_unit='M', t_unit=None, fig_size_scaler=1,\
                      plot_mode='Separate',DP_scaler=1):
    '''Generates plots for each experiment in the kinetic dataset "data". 
     If the stochiometry_list argument is specified, plot_data_MB from auto_vtna 
     is used to also show the mass balance for the data.'''
    if plot_mode=='Together':
        plot_mode='together'
    elif plot_mode=='Scroll':
        plot_mode='scrollable'
    if stochiometry_list != None:
        # Obtain species and stochiometry lists from dict input.
        species = list(stochiometry_list.keys())
        stochiometry = list(float(i) for i in stochiometry_list.values())
        # plot the data with mass balance. 
        auto_vtna.plot_data_MB(data, species, stochiometry, ylim=ylim, t_unit=t_unit, y_unit=y_unit, fig_size_scaler=fig_size_scaler,plot_mode=plot_mode,DP_scaler=DP_scaler)
    else:
        # Plot the data without mass balance shown. 
        auto_vtna.plot_data(data, ylim=ylim, y_unit=y_unit, t_unit=t_unit, fig_size_scaler=fig_size_scaler,plot_mode=plot_mode,DP_scaler=DP_scaler)

def is_float(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

def is_int(value):
    try:
        int(value)
        return True
    except ValueError:
        return False

def run_calculation(data, selected_sheets,selected_columns, selected_result_columns, calculation_settings):
    # Create a VTNA selection dictonary.
    VTNA_selection_dict=auto_vtna.make_VTNA_selection(data)
    # Add the normalised species selected in the GUI with orders of 1 as placeholders.
    VTNA_selection_dict['normalised_species']={}
    for molecule in selected_columns:
        VTNA_selection_dict['normalised_species'][molecule]=1
    # Define the experiments selected in the GUI that will be included in the calculation.
    VTNA_selection_dict['RM']={}
    for i in selected_sheets:
        VTNA_selection_dict['RM'][i]={}
        VTNA_selection_dict['RM'][i]["omissions"]=None
    # Define the output species specified in the GUI.
    VTNA_selection_dict['output_species']=selected_result_columns[0]
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
    if calculation_settings['constraint']=='through origin':
        constraint='origin'
    if calculation_settings['constraint']=='None':
        constraint=None
    # Set up the automatic VTNA calculation. 
    VTNA_auto_calculation=Automatic_VTNA(data,VTNA_selection_dict,order_range=calculation_settings["order_range"],
    resolution=calculation_settings["resolution"],deg=calculation_settings["deg"],fit_metric=calculation_settings["fit_metric"],
    iterations=calculation_settings["iterations"], constraint=constraint,score_interval=calculation_settings["score_interval"],
    fixed_order_species=fixed_order_species)
    return VTNA_auto_calculation

def check_kinetic_data(kinetic_data):
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
    words = sentence.split()
    for i, word in enumerate(words):
        if len(word) > length_cutoff:
            middle_index = len(word) // 2
            words[i] = word[:middle_index] + "- " + word[middle_index:]
    return ' '.join(words)

def apply_best_number_format(x, significant_figures):
    # Check if the initial concentration value can be rounded or if 
    # scientific notation is needed. Apply the right modification. 
    number_legth1=len(str(x))
    number_length2=len(str(round(x,significant_figures)))
    if number_legth1>number_length2:
        return f'{x:.{significant_figures}e}'
    else:
        return round(x,significant_figures)

def generate_initial_concentration_table(kinetic_data,fig_size,conc_unit=None,significant_figures=3):
    # NB: this code was written with help from OpenAI's ChatGTP.
    # Close any existing plots.
    plt.close('all')
    # Obtain the intial concentrations for each experiment in the kinetic dataset. 
    initial_concs_df= auto_vtna.initial_concs(kinetic_data)
    # Format numeric values with scientific notation and specified significant figures
    initial_concs_df_formatted = initial_concs_df.applymap(
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
            if len(initial_concs_df.columns)>6 and len(initial_concs_df.columns)<10:
                text=add_dash_to_long_words(text,11)
            if len(initial_concs_df.columns)>10:
                text=add_dash_to_long_words(text,7)
            text=text.replace(' ','\n')
            table._cells[i].set_text_props(text=text)
        # Set the cells in the first column to blue.
        if i[1]== 0:
            table._cells[i].set_facecolor('#6699CC')
    # Remove the top left corner cell from the table by setting its colour to white.
    table._cells[(0,0)].set_text_props(text='')
    table._cells[(0,0)].set_facecolor('#FFFFFF')
    table._cells[(0,0)].set_linewidth(0)
    plt.tight_layout()
    plt.show()

def plot_concentration_profiles(kinetic_data, selected_species,ylim=None, y_unit='M', t_unit=None, fig_size_scaler=1,DP_scaler=1):
    plt.close('all')
    auto_vtna.plot_data_together(kinetic_data, selected_species,ylim=ylim, y_unit=y_unit, t_unit=t_unit, fig_size_scaler=fig_size_scaler*0.7,DP_scaler=DP_scaler)
    
def main():
    excel_loaded = False  # Flag to track whether an Excel file is loaded
    img = None 
    fixed_orders={}
    enable_plot_button = False  # Flag to track whether the calculation has been run
    calculation_settings = {"order_range": [-1.5, 2.5], "resolution": 10, "deg": 5, "fit_metric": 'RMSE',
                            "iterations": 4, "constraint": 'monotonic', "score_interval": 0.15,
                            "fixed_order_species": None}
    order_vs_overlay_plot_settings = {"y_unit":'M',"size_scaler1":1,"popup_scaler":1,"colour_scaler":1,
                            "contour_resolution":None,"datapoints":False,"interval":False,"zoom":'None',
                            "contour_cbar_max":None,"contour_cbar_min":None, "custom_title":None,"show_legend_popup":True,
                            'show_legend_main':True,'decimals_cbar':3,'annotate':True}
    auto_vtna_results={"best_order_dict":{'':'','\t':'','\r':''},"order_vs_overlay_table":None,"interval_table":None,"plot_image":None}
    data_plotting_settings={"ylim":None,"t_unit":None,'y_unit':'M',"size_scaler":1,'DP_scaler':1,"mode":'together',"significant_figs":3}
    overlay_plot_settings={'y_unit':'M',"size_scaler":1,"tt_scaler":1,"grid":False,"custom_title":None,"check_score":False,\
                           "check_fit":False,"constraint":'monotonic','deg':5,"fit_metric":'RMSE',
                           'score_settings_visibility':False,'DP_scaler':1,'xtick_rotation':0,'legend':True,
                           'legend_position':'Inside'}
    best_order_dict={}
    best_order_dict_rounded={}
    omission_dictionary=None
    omission_range_dictionary={}
    omission_range_dictionary["range_type"]='Percentage'
    omission_mode='Datapoint mode'
    Plot_concentration_profiles_selected_sheets=False
    layout = [
        [sg.Frame("Load Kinetic Data", layout=[[sg.Text("Use Excel File:"),
        sg.InputText(key="file_path", enable_events=True, size=(12, 1)), sg.FileBrowse(),sg.Text("Or:"),sg.Button("Use CSVs",key="Select CSVs")]])],
        [sg.Button("Check Data",key="Check Kinetic Data", disabled=True), sg.Button("Visualise Data",key="Inspect Kinetic Data", disabled=True),
        sg.Button("Crop Data", key="Crop Kinetic Data", disabled=True),sg.Combo(['Range mode', 'Datapoint mode'], default_value=omission_mode,
                            key="crop_mode", enable_events=True, readonly=True,size=(12,1))],
        [sg.Frame("Select Sheets", layout=[
            [sg.Listbox(values=[], select_mode=sg.LISTBOX_SELECT_MODE_MULTIPLE, size=(22, 5), key="sheets")]
        ]),
        sg.Frame("Select Normalised Species", layout=[
            [sg.Listbox(values=[], select_mode=sg.LISTBOX_SELECT_MODE_MULTIPLE, size=(21, 5), key="columns")]
        ])],
        [sg.Frame("Select output species", layout=[
            [sg.Listbox(values=[], select_mode=sg.LISTBOX_SELECT_MODE_SINGLE, size=(20, 1), key="result_columns")]
        ]),
        sg.Frame("VTNA overlay plot", layout=[
            [sg.Button("Generate plot", key="Generate Overlay Plot", disabled=True),
            sg.Button("Plot Settings", key="Overlay Plot Settings", disabled=True)]
        ])],
        [sg.Frame("Automatic VTNA", layout=[
            [sg.Button(" Run ", key='Run Auto VTNA' ,disabled=True),
            sg.Button(" Calculation ⚙", key="calculation_settings", disabled=True),
            sg.Button(" Save \U0001F4BE ", key="Save Results", disabled=not enable_plot_button)],
            [sg.Button("Plot Results", disabled=not enable_plot_button),
            sg.Button("Plot ⚙", key="Plot Settings", disabled=not enable_plot_button),
            sg.Button("Show +Info", key="extra_info", disabled=not enable_plot_button)]
        ]),
        sg.Table(values=[[key, value] for key, value in auto_vtna_results["best_order_dict"].items()],
                headings=['Compound', 'Order'],
                auto_size_columns=False,
                justification='right',
                num_rows=min(5, len(auto_vtna_results["best_order_dict"])),
                col_widths=[8, 5],
                key='-TABLE-')],
        [sg.pin(sg.Text('placeholder',key="Overlay Score",enable_events=True,visible=False,metadata=False,font=11))],
        [sg.pin(sg.Table(values={},headings=['Normalised species','Lower order limit','Upper order limit'],auto_size_columns=False,col_widths=[14,13,13],key='interval_table',visible=False))]
    ]
    window = sg.Window("Automatic VTNA Calculator", layout)
    plt.plot([1,2],[1,2])
    plt.close()
    while True:
        event, values = window.read()
        if event == sg.WIN_CLOSED:
            break
        elif event == "file_path":
            # Update the list of sheets when a file is selected
            file_path = values["file_path"]
            files_calculation=file_path
            if file_path:
                try:
                    kinetic_data = pd.read_excel(file_path, sheet_name=None)
                    # Convert column names to strings.
                    for sheet_name, df in kinetic_data.items():
                        df.columns = df.columns.astype(str)
                    omission_range_dictionary={}
                    omission_range_dictionary["range_type"]='Percentage'
                    data = kinetic_data.copy()
                    df_nrows=[]
                    for i in kinetic_data.keys():
                        df_nrows.append(kinetic_data[i].shape[0])
                    if max(df_nrows)>25:
                        omission_mode='Range mode'
                    else:
                        omission_mode='Datapoint mode'
                    window["crop_mode"].update(value=omission_mode)
                    time_label=list(data.values())[0].columns[0]
                    sheets = list(data.keys())
                    columns = data[sheets[0]].columns[1:].tolist()  # Exclude the first column
                    fixed_order_species = data[sheets[0]].columns[1:].tolist()  # All reaction species
                    window["sheets"].update(values=sheets)
                    window["columns"].update(values=columns)
                    window["result_columns"].update(values=columns)
                    excel_loaded = True  # Set the flag to True
                    window["Run Auto VTNA"].update(disabled=False)  # Enable the "Run Automatic VTNA" button
                    window["calculation_settings"].update(disabled=False)  # Enable the "Calculation Settings" button
                    window["Check Kinetic Data"].update(disabled=False)
                    window["Inspect Kinetic Data"].update(disabled=False)
                    window["Generate Overlay Plot"].update(disabled=False)
                    window["Overlay Plot Settings"].update(disabled=False)
                    window["Crop Kinetic Data"].update(disabled=False)
                except PermissionError as e:
                    sg.popup_error(f"Error opening file (make sure it is not open): {e}")
                    continue
        elif event == "crop_mode":
            omission_mode=values["crop_mode"]
        elif event == "Save Results":
            layout_save = [
                [sg.Text("Filename:"), sg.InputText(key="filename", size=(5, 1))],
                [sg.Text("File type:"), sg.Combo(['csv', 'xlsx'], default_value='xlsx', key="filetype", readonly=True)],
                [sg.Button("OK"), sg.Button("Cancel")]
            ]
            window_save = sg.Window("Save file settings", layout_save)

            while True:
                event_save, values_save = window_save.read()
                if event_save == sg.WINDOW_CLOSED or event_save == 'OK' or event_save == 'Cancel':
                    if window_save['filename'] is None or window_save['filename'] == '':
                        raise ValueError("Define a filename to continue.")
                        break   
                    filetype = values_save["filetype"]
                    filename = values_save['filename']
                    files_calculation_nested = [['', files_calculation]]

                    if type(files_calculation) != str:
                        files_calculation_nested = [['', i] for i in files_calculation]

                    files_calculation_nested.insert(0, ["Kinetic data is loaded from:"])
                    fixed_order_species_name=VTNA_auto_calculation.fixed_order_species
                    if fixed_order_species_name==None:
                        fixed_order_species_name='None'
                    constraint = VTNA_auto_calculation.constraint if VTNA_auto_calculation.constraint not in [None,''] else 'None'
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
                        ["  - Fixed order species:", fixed_order_species_name],
                        ["VTNA SELECTION DICTIONARY:"]
                    ]
                    selection_list = list(VTNA_auto_calculation.VTNA_selection.items())
                    s0_2=[['RM']]
                    for key, value in selection_list[0][-1].items():
                        RM_omissions=[key]
                        if list(list(value.items())[0])[1]!=None:
                            RM_omissions=RM_omissions+[list(list(value.items())[0])[0]]+list(list(value.items())[0])[1]
                        else:
                            RM_omissions=RM_omissions+[list(list(value.items())[0])[0]]+['None']
                        s0_2.append(RM_omissions)
                    s1_2 = [['normalised_species']]
                    s1_2.extend(["  - "+key, value] for key, value in selection_list[2][1].items())
                    selection_list[2] = s1_2
                    selection_list=s0_2+[selection_list[1]]+s1_2
                    page1 = files_calculation_nested + page1 + selection_list
                    page1_df = pd.DataFrame(page1)
                    if filetype == 'csv':
                        try:
                            with open(f'{filename}.csv', 'w', newline='') as csvfile:
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
                            with pd.ExcelWriter(f'{filename}.xlsx', engine='xlsxwriter') as writer:
                                page1_df.to_excel(writer, index=False, header=False, sheet_name='Calculation settings')

                            with pd.ExcelWriter(f'{filename}.xlsx', engine='openpyxl', mode='a') as writer:
                                VTNA_auto_calculation.interval.to_excel(writer, index=False, header=True, sheet_name=f'{str(float(VTNA_auto_calculation.score_interval)*100)}% order intervals.')

                            with pd.ExcelWriter(f'{filename}.xlsx', engine='openpyxl', mode='a') as writer:
                                VTNA_auto_calculation.results.to_excel(writer, index=False, header=True, sheet_name='Overlay score vs. order matrix')

                            if img is not None:
                                writer.book.create_sheet('Overlay score vs. order plot')
                                writer.book.save(f'{filename}.xlsx')
                                sheet = writer.sheets['Overlay score vs. order plot']
                                sheet.add_image(img, 'A1')
                                writer.book.save(f'{filename}.xlsx')
                        except Exception as e:
                            sg.popup_error(f"Filename permission error (close excel file with same name): {e}")
                            continue
                    window_save.close()
                    break
        elif event == "Run Auto VTNA":
            # Get selected sheet, columns, and result columns
            selected_sheet = values["sheets"]
            selected_columns = values["columns"]
            selected_result_columns = values["result_columns"]
            if selected_sheet and selected_columns and selected_result_columns:
                count=len(selected_columns)
                if type(calculation_settings["fixed_order_species"])==dict:
                    for i in selected_columns:
                        if i in list(calculation_settings["fixed_order_species"].keys()):
                            count+=-1
                if count>=1:
                        VTNA_auto_calculation = run_calculation_in_thread(window, data, selected_sheet, selected_columns, selected_result_columns, calculation_settings)
                        best_order_dict={}
                        for i in range(len(VTNA_auto_calculation.best_orders)):
                            best_order_dict[f"{VTNA_auto_calculation.results.columns[i].replace('order in ','')}"]=VTNA_auto_calculation.best_orders[i]
                        orders_last_species=VTNA_auto_calculation.results.iloc[:,-2].to_numpy()
                        resolution_orders=str(format(orders_last_species[-1]-orders_last_species[-2],".12f"))
                        best_order_dict_rounded={}
                        for i in range(len(resolution_orders)):
                            if resolution_orders[i]!='.' and resolution_orders[i]!='0':
                                break
                        for species,best_order in best_order_dict.items():
                            best_rounded=f'{round(best_order, i-1):.{int(i-1)}f}'
                            best_order_dict_rounded[species]=best_rounded
                        auto_vtna_results["best_order_dict"]=best_order_dict
                        auto_vtna_results["order_vs_overlay_table"]=VTNA_auto_calculation.results
                        auto_vtna_results["interval_table"]=VTNA_auto_calculation.interval
                        if VTNA_auto_calculation:
                            sg.popup("Calculation complete.",auto_close=True, auto_close_duration=2)
                            window['-TABLE-'].update(values=[[key, value] for key, value in best_order_dict_rounded.items()])
                            enable_plot_button = True  # Enable the "Plot Results" button
                            # Update the button state
                            window["Plot Results"].update(disabled=not enable_plot_button)
                            window["Save Results"].update(disabled=not enable_plot_button)
                            window["Plot Settings"].update(disabled=not enable_plot_button)
                            window["extra_info"].update(disabled=False)
                            if calculation_settings['fit_metric']=='R2':
                                best_OS=round(max(VTNA_auto_calculation.results[f"{VTNA_auto_calculation.fit_metric} overlay score"]),5)
                            else:
                                best_OS=round(min(VTNA_auto_calculation.results[f"{VTNA_auto_calculation.fit_metric} overlay score"]),6)  
                            overlay_score_string=f"Overlay score ({calculation_settings['fit_metric']}) at the optimal order values: {str(best_OS)} \nIntervals of order values less than {round(calculation_settings['score_interval']*100)}% away from optimum:"   
                            df_interval=VTNA_auto_calculation.interval.round(decimals=5)
                            table_values=df_interval.values.tolist()
                            for i,value in enumerate(table_values):
                                table_values[i]=tuple(value)
                            window['interval_table'].update(values=table_values)
                            window["Overlay Score"].update(overlay_score_string)
                            window['interval_table'].update(num_rows=len(table_values))   
                else:
                    sg.popup_error("The number of fixed order species can't be the same as or higher than the number of normalised reaction species.")
            else:
                sg.popup_error("Please select sheets, reaction species, and output reaction species to run automatic VTNA.")
        elif event == 'extra_info':
            # Access the me
            element_table = window['interval_table']
            visible = element_table.metadata = not element_table.metadata
            window['interval_table'].update(visible=visible)
            window["Overlay Score"].update(visible=visible)
            if visible:
                window['extra_info'].update('Hide + Info')
            else:
                window['extra_info'].update('Show +Info')

        elif event == "Plot Settings":
            layout_settings_OvsO_plot = [
                [sg.Frame("General settings",layout=[
                    [sg.Text("Show legend in main plot:"),
                    sg.Combo(['True','False'], default_value=str(order_vs_overlay_plot_settings['show_legend_main']),key='show_legend_main', enable_events=True, readonly=True,size=(6,1))],
                    [sg.Text("Show error interval datapoints:"),
                    sg.Combo(['True','False'], default_value=str(order_vs_overlay_plot_settings['interval']),key="interval", enable_events=True, readonly=True,size=(6,1))],
                    [sg.Text("Overlay score range:"),sg.InputText(key="contour_cbar_min", size=(6, 1),default_text='Default' if order_vs_overlay_plot_settings["contour_cbar_min"] is None else str(order_vs_overlay_plot_settings["contour_cbar_min"])),
                    sg.Text(" to "),sg.InputText(key="contour_cbar_max", size=(6, 1),default_text='Default' if order_vs_overlay_plot_settings["contour_cbar_max"] is None else str(order_vs_overlay_plot_settings["contour_cbar_max"]))],
                    [sg.Text("Custom figure title:"),
                    sg.InputText(key="custom_title", size=(19, 1),default_text=str(order_vs_overlay_plot_settings["custom_title"]) if order_vs_overlay_plot_settings["custom_title"]!=None else '')],
                    [sg.Text("Main figure size scaler:    "),
                    sg.InputText(key="main_size_scaler", size=(5, 1),default_text=float(order_vs_overlay_plot_settings["size_scaler1"]))]])],
                [sg.Frame("Settings for contour plot",layout=[
                    [sg.Text("Show datapoints:            "),
                     sg.Combo(['True','False'], default_value=str(order_vs_overlay_plot_settings['datapoints']),key="datapoints", enable_events=True, readonly=True,size=(6,1))],
                    [sg.Text("Show optimum annotation:"),
                     sg.Combo(['True','False'], default_value=str(order_vs_overlay_plot_settings['annotate']),key="annotate", enable_events=True, readonly=True,size=(6,1))],
                    [sg.Text("Color bar max scaler:      "),
                    sg.InputText(key="colour_scaler", size=(6, 1),default_text=str(order_vs_overlay_plot_settings["colour_scaler"]))],
                    [sg.Text("Colour bar resolution:      "),
                    sg.InputText(key="contour_resolution", size=(6, 1),default_text='Default' if order_vs_overlay_plot_settings["contour_resolution"] is None else str(order_vs_overlay_plot_settings["contour_resolution"]))],
                    [sg.Text("Domain around optimum (zoom):"),
                     sg.InputText(key="zoom", size=(5, 1),default_text=str(order_vs_overlay_plot_settings["zoom"]))],
                    [sg.Text("Number of decimals for cbar:   "),
                     sg.InputText(key="decimals_cbar", size=(5, 1),default_text=str(order_vs_overlay_plot_settings["decimals_cbar"]))]])],
                [sg.Frame("Settings for pop-up overlay plots",layout=[
                    [sg.Text("Show legend:            "),
                    sg.Combo(['True','False'], default_value=str(order_vs_overlay_plot_settings['show_legend_popup']),key="show_legend_popup", enable_events=True, readonly=True,size=(6,1))],
                    [sg.Text("Pop-up figure size scaler:"),
                    sg.InputText(key="popup_size_scaler", size=(6, 1),default_text=str(order_vs_overlay_plot_settings["popup_scaler"]),enable_events=True)],
                    [sg.Text("Concentration unit:          "),
                    sg.InputText(key="y_unit", size=(6, 1),default_text=str(order_vs_overlay_plot_settings["y_unit"]))]])],
                [sg.Button("OK"), sg.Button("Cancel")]]
                    
            window_OvsO_plot_settings = sg.Window("Order Versus Overlay Plot Settings", layout_settings_OvsO_plot)
            while True:
                event_OvsO_plot_settings, values_OvsO_plot_settings = window_OvsO_plot_settings.read()
                if event_OvsO_plot_settings == sg.WIN_CLOSED or event_OvsO_plot_settings == "Cancel":
                    break
                elif event_OvsO_plot_settings == "OK":
                    try:
                        order_vs_overlay_plot_settings["y_unit"] = str(values_OvsO_plot_settings["y_unit"])
                        main_size_scaler = float(values_OvsO_plot_settings["main_size_scaler"])
                        popup_scaler = float(values_OvsO_plot_settings["popup_size_scaler"])
                        colour_scaler = float(values_OvsO_plot_settings["colour_scaler"])
                        if order_vs_overlay_plot_settings["contour_resolution"] not in [None,'Default']:
                            order_vs_overlay_plot_settings["contour_resolution"]=float(values_OvsO_plot_settings["contour_resolution"])
                        else:
                            order_vs_overlay_plot_settings["contour_resolution"]=None
                        datapoints = values_OvsO_plot_settings["datapoints"]
                        annotate = values_OvsO_plot_settings["annotate"]
                        interval = values_OvsO_plot_settings["interval"]
                        show_legend_popup = values_OvsO_plot_settings["show_legend_popup"]
                        show_legend_main = values_OvsO_plot_settings["show_legend_main"]
                        zoom_range = str(values_OvsO_plot_settings["zoom"])
                        decimals_cbar = values_OvsO_plot_settings["decimals_cbar"]
                        order_vs_overlay_plot_settings["decimals_cbar"]=decimals_cbar
                        if values_OvsO_plot_settings["contour_cbar_max"] not in [None,'Default']:
                            order_vs_overlay_plot_settings["contour_cbar_max"] = float(values_OvsO_plot_settings["contour_cbar_max"])
                        else:
                            order_vs_overlay_plot_settings["contour_cbar_max"] = None
                        if values_OvsO_plot_settings["contour_cbar_min"] not in [None,'Default']:
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
                        # Check size scaler value
                        if not (0.1 < main_size_scaler <= 3):
                            raise ValueError("Invalid figure size scaler. Must be between 0.1 and 3.")
                        if not (0.1 < popup_scaler <= 3):
                            raise ValueError("Invalid popup figure size scaler. Must be between 0.1 and 3.")
                        #if not (0 < colour_scaler <= 1):
                            #raise ValueError("Invalid contour plot colour bar scaler. Must be between 0 and 1")
                        if zoom_range!='None':
                            if not (0 < float(zoom_range) <= float(calculation_settings["order_range"][-1]-calculation_settings["order_range"][0])+2):
                                raise ValueError("The zoom range defines the breadth of order values shown in the\
 order versus overlay plot, and can't be below 0 and should't be much higher than the order range of the calculation")
                        order_vs_overlay_plot_settings["size_scaler1"] = main_size_scaler
                        order_vs_overlay_plot_settings["popup_scaler"] = popup_scaler
                        order_vs_overlay_plot_settings["colour_scaler"] = colour_scaler
                        order_vs_overlay_plot_settings["custom_title"]= str(values_OvsO_plot_settings["custom_title"])
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
            # Open Settings window
            if excel_loaded:
                layout_settings = [
                    [sg.Text("Overlay Plot Settings", font=("helvetica", 12, "bold"))],
                    [sg.Text("Concentration unit:       "),
                    sg.InputText(key="y_unit", size=(5, 1), default_text=str(overlay_plot_settings["y_unit"]))],
                    [sg.Text("tT axis label scaler:      "),
                    sg.InputText(key="tt_scaler", size=(5, 1), default_text=str(overlay_plot_settings["tt_scaler"]))],
                    [sg.Text("Datapoint size scaler:   "),
                    sg.InputText(key="DP_scaler", size=(5, 1), default_text=str(overlay_plot_settings["DP_scaler"]))],
                    [sg.Text("x-axis tick label rotation:"),
                    sg.InputText(key="xtick_rotation", size=(5, 1), default_text=str(overlay_plot_settings["xtick_rotation"]))],
                    [sg.Text("Figure size scaler:        "),
                    sg.InputText(key="size_scaler", size=(5, 1), default_text=str(overlay_plot_settings["size_scaler"]))],
                    [sg.Text("Show gridlines:            "),
                    sg.Combo(['True', 'False'],size=(5, 1), default_value=str(overlay_plot_settings['grid']), key="grid", enable_events=True, readonly=True)],
                    [sg.Text("Show legend:              "),
                    sg.Combo(['True', 'False'],size=(5, 1), default_value=str(overlay_plot_settings['legend']), key="legend", enable_events=True, readonly=True)],
                    [sg.Text("Legend position:          "),
                    sg.Combo(['Inside','Outside'],size=(5, 1), default_value=str(overlay_plot_settings['legend_position']), key="legend_position", enable_events=True, readonly=True)],
                    [sg.Text("Custom title:"),
                    sg.InputText(key="custom_title", size=(22, 1),default_text=str(overlay_plot_settings["custom_title"]) if overlay_plot_settings["custom_title"]!=None else 'Default')],
                    [sg.Checkbox("Calculate overlay score", key='check_score',default=overlay_plot_settings["check_score"],enable_events=True),
                     sg.Checkbox("Show fit function", key='check_fit',default=overlay_plot_settings["check_fit"],enable_events=True)],
                    [sg.pin(sg.Text("Constraint:",enable_events=True,visible=overlay_plot_settings['score_settings_visibility'],key="constraint_text")),sg.pin(sg.Combo(['monotonic', 'None', 'through origin'],
                            default_value=overlay_plot_settings["constraint"],
                            key="constraint", enable_events=True, readonly=True,visible=overlay_plot_settings['score_settings_visibility']))],
                    [sg.pin(sg.Text("Degree:",enable_events=True,visible=overlay_plot_settings['score_settings_visibility'],key='deg_text')),sg.pin(sg.Combo([1, 2, 3, 4, 5, 6, 7], default_value=overlay_plot_settings["deg"],
                            key="deg",enable_events=True, readonly=True,visible=overlay_plot_settings['score_settings_visibility']))],
                    [sg.pin(sg.Text("Goodness-of-Fit Measurement:",enable_events=True,visible=overlay_plot_settings['score_settings_visibility'],key='GOF_text')),sg.pin(sg.Combo(['SE', 'RMSE', 'R2', 'Variance'], 
                            default_value=overlay_plot_settings["fit_metric"],
                            key="fit_metric", enable_events=True, readonly=True,visible=overlay_plot_settings['score_settings_visibility']))],
                    [sg.Button("OK"), sg.Button("Cancel")]
                ]
                window_overlay_settings = sg.Window("Overlay Plot Settings", layout_settings,resizable=True)
                while True:
                    event_overlay_settings, values_overlay_settings = window_overlay_settings.read()
                    if event_overlay_settings == sg.WIN_CLOSED or event_overlay_settings == "Cancel":
                        break
                    elif event_overlay_settings in ['check_score','check_fit']:
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
                    elif event_overlay_settings == "OK":
                        try:
                            y_unit = values_overlay_settings["y_unit"]
                            tt_scaler = float(values_overlay_settings["tt_scaler"])
                            DP_scaler = float(values_overlay_settings["DP_scaler"])
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
                            deg=values_overlay_settings["deg"]
                            constraint=values_overlay_settings["constraint"]
                            overlay_plot_settings["fit_metric"]=values_overlay_settings["fit_metric"]
                            overlay_plot_settings["check_score"]=values_overlay_settings["check_score"]
                            overlay_plot_settings["check_fit"]=values_overlay_settings["check_fit"]
                            if grid=='True':
                                overlay_plot_settings["grid"]=True
                            else:
                                overlay_plot_settings["grid"]=False
                            if legend_bol=='True':
                                overlay_plot_settings["legend"]=True
                            else:
                                overlay_plot_settings["legend"]=False
                            # Check size scaler value
                            if deg not in [1, 2, 3, 4, 5, 6, 7,8,9,10,11,12,13]:
                                raise ValueError("Invalid degree. Must be an integer between 1 and 13.")
                            if constraint=='through origin':
                                if int(deg)!=1:
                                    raise ValueError("Set fitting degree to 1 for linear line through origin.")
                            if not (0.1 < size_scaler <= 3):
                                raise ValueError("Invalid figure size scaler. Must be between 0.1 and 3.")
                            if not (0.1 < tt_scaler <= 3):
                                raise ValueError("Invalid figure transformed time axis label scaler. Must be between 0.1 and 3.")
                            if not (0.00 <= DP_scaler <= 5):
                                raise ValueError("Invalid datapoint size scaler. Must be between 0 and 5.")
                            if not (0.00 <= xtick_rotation <= 180):
                                raise ValueError("Invalid x-axis tick label rotation value. Must be between 0 and 180.")
                            overlay_plot_settings["deg"]=deg
                            overlay_plot_settings["constraint"]=constraint
                            overlay_plot_settings["xtick_rotation"]=xtick_rotation
                        except ValueError as e:
                                sg.popup_error(f"Invalid input: {e}")
                                continue
                        window_overlay_settings.close()
                window_overlay_settings.close()
        elif event == "Plot Results":
            try:
                no_dimensions=-len(fixed_orders.values())+len(selected_columns)
                if no_dimensions>2:
                    raise ValueError("Can't plot results for more than two normalised reaction species to investigate.")
                plt.close('all')
                # Run the plotting based on the calculation results
                zoom_range=order_vs_overlay_plot_settings["zoom"]
                if zoom_range=='None':
                    zoom_range=False
                custom_title = order_vs_overlay_plot_settings["custom_title"]
                if custom_title!=None and custom_title!='Default':
                    if len(custom_title)==0 or custom_title=='None':
                        custom_title=None
                img=VTNA_auto_calculation.plot_orders_vs_overlay(y_unit=order_vs_overlay_plot_settings["y_unit"],size_scaler=order_vs_overlay_plot_settings["size_scaler1"],
                size_scaler2=order_vs_overlay_plot_settings["popup_scaler"],color_scaler=order_vs_overlay_plot_settings["colour_scaler"],
                fixed_cbar_resolution=order_vs_overlay_plot_settings["contour_resolution"],points=order_vs_overlay_plot_settings["datapoints"],
                interval=order_vs_overlay_plot_settings["interval"],zoom=zoom_range,overlay_score_range_max=order_vs_overlay_plot_settings["contour_cbar_max"],
                overlay_score_range_min=order_vs_overlay_plot_settings["contour_cbar_min"],title=custom_title,show_legend_popup=order_vs_overlay_plot_settings["show_legend_popup"],
                show_legend_main=order_vs_overlay_plot_settings["show_legend_main"],decimals_in_colorbar=order_vs_overlay_plot_settings["decimals_cbar"],
                annotate_optimum=order_vs_overlay_plot_settings["annotate"])
            except ValueError as e:
                sg.popup_error(f"Invalid input: {e}")
                continue
        elif event == "Generate Overlay Plot":
            plt.close('all')
            selected_sheet = values["sheets"]
            selected_columns = values["columns"]
            selected_result_columns = values["result_columns"]
            tt_scaler=overlay_plot_settings['tt_scaler']
            DP_scaler=float(overlay_plot_settings["DP_scaler"])
            xtick_rotation=float(overlay_plot_settings["xtick_rotation"])
            y_unit=overlay_plot_settings['y_unit']
            size_scaler=overlay_plot_settings['size_scaler']
            grid=overlay_plot_settings['grid']
            title=None if overlay_plot_settings['custom_title']=='Default' else overlay_plot_settings['custom_title']
            constraint=overlay_plot_settings['constraint']
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
            if selected_sheet and selected_columns and selected_result_columns:
                # Get selected sheet, columns, and result columns
                selected_sheet = values["sheets"]
                selected_columns = values["columns"]
                selected_result_columns = values["result_columns"]
                # Open a pop-up window to input order values for selected fixed order species
                if selected_columns:
                    layout_order_values = [
                        [sg.Text(f"Enter order values for selected species:")],
                    ]
                    for species in selected_columns:
                        layout_order_values.append(
                            [sg.Text(f"{species}"), sg.InputText(key=f"order_value_{species}", size=(5, 1),enable_events=True)])
                    buttons=[sg.Button("OK"), sg.Button("Cancel")]
                    if len(best_order_dict)!=0:
                        if [i in selected_columns for i in list(best_order_dict.keys())]:
                            buttons.append(sg.Button("Auto-fill calculated values.",key='auto'))
                    layout_order_values.append(buttons)
                    window_order_values = sg.Window("Enter Order Values", layout_order_values,resizable=True,finalize=True)
                    while True:
                        input_valid = False
                        event_order_values, values_order_values = window_order_values.read()
                        if event_order_values=='auto':
                            for species,order in best_order_dict_rounded.items():
                                window_order_values[f"order_value_{species}"].update(value=order)
                        elif event_order_values == sg.WIN_CLOSED or event_order_values == "Cancel":
                            # User closed the window or pressed "Cancel"
                            plt.close('all')
                            window_order_values.close()
                            break
                        if event_order_values=="OK":
                            order_values = {}
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
                                overlay_plot(data,selected_sheet,order_values,selected_result_columns[0],tt_scaler=tt_scaler,\
                                         grid=grid,y_unit=y_unit,size_scaler=size_scaler,title=title, fit_metric=fit_metric,
                                         deg=deg,constraint=constraint,show_overlay_score=show_overlay_score,\
                                         show_fit_function=show_fit_function,DP_scaler=DP_scaler,xtick_rotation=xtick_rotation,\
                                         show_legend=legend_bol,legend_outside=legend_outside)
            else:
                sg.popup_error("Please select sheets, reaction species, and output reaction species before generating overlay plot.")
        elif event == "calculation_settings":
            # Open Calculation Settings window
            if excel_loaded:
                fixed_orders={}
                layout_fixed_order=[]
                label_width = max(len(str(compound)) for compound in fixed_order_species)
                fixed_orders_start=calculation_settings["fixed_order_species"]
                if type(fixed_orders_start)!=dict:
                    fixed_orders_start={}
                for compound in fixed_order_species:
                    if compound in fixed_orders_start.keys():
                            layout_fixed_order.append([
                                sg.Text(compound, size=(label_width, 1)),
                                sg.Checkbox("", key=f"Check_{compound}", enable_events=True,default=True),
                                sg.InputText(key=f"Value_{compound}", visible=True,size=(5,1),
                                default_text=str(fixed_orders_start[compound]))])
                    else:
                        layout_fixed_order.append([
                            sg.Text(compound, size=(label_width, 1)),
                            sg.Checkbox("", key=f"Check_{compound}", enable_events=True),
                            sg.InputText(key=f"Value_{compound}", visible=False,size=(5,1),
                            )])
                layout_settings = [
                [sg.Text("Calculation Settings", font=("helvetica", 12, "bold"))],

                [sg.Frame("Fit Settings", layout=[
                    [sg.Text("Constraint:"),
                    sg.Combo(['monotonic', 'None', 'through origin'],
                            default_value=calculation_settings["constraint"],
                            key="constraint", enable_events=True, readonly=True,size=(11,1)),
                    sg.Text("Degree:"),
                    sg.Combo([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], default_value=calculation_settings["deg"], key="deg")],
                    [sg.Text("Goodness-of-fit metric:"),
                    sg.Combo(['SE', 'RMSE', 'R2', 'Variance'], default_value=calculation_settings["fit_metric"],
                            key="fit_metric", enable_events=True, readonly=True,size=(8,1))]
                ])],
                [sg.Frame("Order exploration settings", layout=[
                    [sg.Text("Iterations:"),
                    sg.Combo([1, 2, 3, 4, 5, 6, 7], default_value=calculation_settings["iterations"], key="iterations"),
                    sg.Text("Order grid density:"),
                    sg.InputText(key="resolution", size=(3, 1),
                                default_text=str(calculation_settings["resolution"]))]
                ])],

                [sg.Frame("Order range to explore", layout=[
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
                [sg.Button("OK"), sg.Button("Cancel")]]

                window_settings = sg.Window("Calculation Settings", layout_settings)

                while True:
                    event_settings, values_settings = window_settings.read()

                    if event_settings == sg.WIN_CLOSED or event_settings == "Cancel":
                        break
                    elif event_settings.startswith("Check_"):
                        window_settings[event_settings.replace('Check_','Value_')].update(visible=values_settings[event_settings])
                        if values_settings[event_settings]==False:
                            window_settings[event_settings.replace('Check_','Value_')].update(value='')
                    elif event_settings == "OK":
                        # Validate input values
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
                            if constraint not in ['monotonic', 'None', 'through origin']:
                                raise ValueError("Invalid constraint. Must be 'monotonic', 'None', or 'through origin'.")

                            # Check score interval
                            if not (0 <= score_interval <= 2):
                                raise ValueError("Invalid score interval. Must be between 0 and 2.")
                            # Check resolution
                            if not (5 <= resolution <= 30):
                                raise ValueError("Invalid score interval. Must be between 5 and 30.")
                            # Check degree
                            if deg not in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10,11, 12, 13]:
                                raise ValueError("Invalid degree. Must be an integer between 1 and 13.")
                            if constraint=='through origin':
                                if int(deg)!=1:
                                    raise ValueError("Set fitting degree to 1 for linear line through origin.")
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
                        calculation_settings["resolution"]= round(resolution)
                        calculation_settings["fixed_order_species"]=fixed_orders
                        if len(fixed_orders)==0:
                            calculation_settings["fixed_order_species"]=None
                        break
                window_settings.close()
        elif event == "Select CSVs":
            # Define counter to know which browse and inputtext elements to make visible. 
            count=1
            sg.theme("DefaultNoMoreNagging")
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
                [sg.Button("Add File", key="add_file"), sg.Button("OK")]
            ]
            # Create the window.
            window_csv = sg.Window("CSV Loader", layout)
            pattern = re.compile(r'/([^/]+)\.csv$')
            while True:
                event, values = window_csv.read()
                # If the window is closed, collect the filenames selected as a list.
                if event == sg.WIN_CLOSED or event == "OK":
                    file_paths_full=[]
                    file_names=[]
                    for key,value in values.items():
                        if value==None:
                            break
                        elif "file_name" in key and len(value)!=0:
                            file_paths = value.split(';')
                            file_paths_full=file_paths_full+file_paths
                            index=''.join(char for char in key if char.isdigit())
                            name=values[f"nickname_{int(index)}"]
                            for i in file_paths:
                                if len(name)==0 or len(file_paths)!=1:
                                    name_search = pattern.search(i)
                                    name=name_search.group(1)
                                    file_names.append(name)
                                else:
                                    file_names.append(name)
                    file_paths_no_duplicates,file_names_no_duplicates=[],[]
                    for count,(i,j) in enumerate(zip(file_paths_full,file_names)):
                        if i not in file_paths_full[:count]:
                            file_paths_no_duplicates.append(i)
                            file_names_no_duplicates.append(j)
                    kinetic_data={}
                    for (i,j) in zip(file_paths_no_duplicates,file_names_no_duplicates):
                        dataset=pd.read_csv(f"{i}")
                        # Convert column titles to strings. 
                        dataset.columns = dataset.columns.astype(str)
                        kinetic_data[j]=dataset
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
                    time_label=list(data.values())[0].columns[0]
                    sheets = list(data.keys())
                    columns = data[sheets[0]].columns[1:].tolist()  # Exclude the first column
                    fixed_order_species = data[sheets[0]].columns[1:].tolist()  # All reaction species
                    window["sheets"].update(values=sheets)
                    window["columns"].update(values=columns)
                    window["result_columns"].update(values=columns)
                    excel_loaded = True  # Set the flag to True
                    window["Run Auto VTNA"].update(disabled=False)  # Enable the "Run Automatic VTNA" button
                    window["calculation_settings"].update(disabled=False)  # Enable the "Calculation Settings" button
                    window["Check Kinetic Data"].update(disabled=False)
                    window["Inspect Kinetic Data"].update(disabled=False)
                    window["Generate Overlay Plot"].update(disabled=False)
                    window["Overlay Plot Settings"].update(disabled=False)
                    window["Crop Kinetic Data"].update(disabled=False)
                    files_calculation=file_paths_no_duplicates
                    break
                # If a file is browsed, make the next InputText and FilesBrowser visible. 
                elif 'file_name' in event and event[-1]==str(count):
                    file_paths = values[event].split(';')
                    if len(file_paths)>1:
                        window_csv[f"nickname_{count}"].update(visible=False)
                        window_csv[f"text_{count}"].update(visible=False)
                    count+=1
                    window_csv[f"file_name_{count}"].update(visible=True)
                    window_csv[f"browse_{count}"].update(visible=True)
                    window_csv[f"text_{count}"].update(visible=True)
                    window_csv[f"nickname_{count}"].update(visible=True)
            window_csv.close()
        elif event == "Check Kinetic Data":
            # Placeholder for checking kinetic data
            kinetic_data_result = check_kinetic_data(data)
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
                        # Input boxes for removing data points by applying ranges. 
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
                remove_data_selective=[sg.Frame("Apply cropping separately for each experiment data:",layout=
                layout_frame_selective,visible=True,key='selective_frame')]
                if len(omission_range_dictionary)>1:
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
                rows_range=[remove_data_selective]+[remove_data_total]
                rows_range.append([sg.Text("Relative or absolute times"),
                                    sg.Combo(['Absolute','Percentage'], default_value=omission_range_dictionary["range_type"],
                                            key="range_type", enable_events=True, readonly=True)])
                layout_range = rows_range + [[[sg.Button('OK')],sg.Button("Reset")]]
                window_range = sg.Window('Crop Data Range', layout_range)
                
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
                            print("info",range_info)
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
                        # 
                        if True in ['dense' in i for i in list(values_range_cleaned.keys())]:
                            # Remove datapoints from data_cropping at regular intervals inputed:
                            selection_dictionary_lower_density=auto_vtna.make_VTNA_selection(kinetic_data).copy()
                            sparser_total=1
                            if 'total_less_dense' in list(values_range_cleaned.keys()):
                                sparser_total=abs(int(float(values_range_cleaned['total_less_dense'])))
                            for i in selection_dictionary_lower_density['RM'].keys():
                                sparser=sparser_total
                                if f'{i}_less_dense' in values_range_cleaned.keys():
                                    sparser=abs(int(float(values_range_cleaned[f'{i}_less_dense'])))
                                length=len(kinetic_data[i].iloc[:, 0])
                                all_DP=list(range(0,length))
                                to_keep=all_DP[0::sparser]
                                if to_keep[-1]!=all_DP[-1]:
                                    to_keep.append(all_DP[-1])
                                to_remove=list(set(all_DP) - set(to_keep))
                                # Don't define omissions datapoints if none are to be removed.
                                if len(to_remove)==0:
                                    continue
                                selection_dictionary_lower_density['RM'][i]['omissions']=to_remove
                            data_cropping=VTNA_omissions(data_cropping,selection_dictionary_lower_density)
                                            # make all values in range_info numerical.
                        range_info={key: [float(value) for value in values] for key, values in range_info.items()}
                        selection_dictionary_ranges=auto_vtna.make_VTNA_selection(kinetic_data).copy()
                        range_list=['range']
                        for i in kinetic_data.keys():
                            if i in list(range_info.keys()):
                                selection_dictionary_ranges['RM'][i]['omissions']=range_list+range_info[i][:2]
                            else:
                                if 'total' in list(range_info.keys()):
                                    selection_dictionary_ranges['RM'][i]['omissions']=range_list+range_info['total'][:2]
                        data_cropping=VTNA_omissions(data_cropping,selection_dictionary_ranges,range_mode=values_range["range_type"])
                        if 4 in [len(i) for i in range_info.values()]:
                            selection_dictionary_ranges_2=auto_vtna.make_VTNA_selection(kinetic_data).copy()
                            for i in kinetic_data.keys():
                                if i in list(range_info.keys()):
                                    if len(range_info[i])==4:
                                        selection_dictionary_ranges_2['RM'][i]['omissions']=range_list+range_info[i][2:]
                                else:
                                    if 'total' in list(range_info.keys()):
                                        if len(range_info['total'])==4:
                                            selection_dictionary_ranges_2['RM'][i]['omissions']=range_list+range_info['total'][2:]
                            data_cropping=VTNA_omissions(data_cropping,selection_dictionary_ranges_2)
                        data=data_cropping
                        omission_range_dictionary=values_range
                        window_range.close()
                    if event_range == sg.WINDOW_CLOSED:
                        window_range.close()
                        break
            else:
                rows = []
                if omission_dictionary!=None:
                    for key, df in kinetic_data.items():
                        row_title = [sg.Text(f'Row {key}')]
                        checkboxes = [sg.Checkbox(str(i), key=f'CHECK_{key}__{i}',default=omission_dictionary[f'CHECK_{key}__{i}']) for i in range(1, len(df[time_label]) + 1)]
                        rows.append(row_title + checkboxes)
                else:
                    for key, df in kinetic_data.items():
                        row_title = [sg.Text(f'Row {key}')]
                        checkboxes = [sg.Checkbox(str(i), key=f'CHECK_{key}__{i}') for i in range(1, len(df[time_label]) + 1)]
                        rows.append(row_title + checkboxes)
                layout = rows + [[sg.Button("Reset")],[sg.Button('OK')]]
                window_omissions = sg.Window('Select rows to reversibly remove from kinetic data.', layout)
                while True:
                    event_omissions, values_omissions = window_omissions.read()
                    if event_omissions=="Reset":
                        for key, df in kinetic_data.items():
                            for i in range(1, len(df[time_label])):
                                window_omissions[f'CHECK_{key}__{i}'].update(value=False)
                    omission_dictionary=values_omissions
                    selection_dictionary=auto_vtna.make_VTNA_selection(kinetic_data)
                    dic=selection_dictionary.copy()
                    for i in dic['RM'].keys():
                        rm_list=[]
                        for j in values_omissions.keys():
                            if i in j and values_omissions[j]:
                                rm_list.append(int(j[-3:].replace('_',''))-1)
                        if len(rm_list)>0:
                            dic['RM'][i]['omissions']=rm_list
                    data=VTNA_omissions(kinetic_data,dic)
                    if event_omissions == sg.WINDOW_CLOSED or event_omissions == 'OK':
                        window_omissions.close()
                        break
        elif event == "Inspect Kinetic Data":
            plt.close('all')
            # Open a new window for inspecting kinetic data
            layout_inspect_kinetic_data = [
                [sg.Button("Generate Initial Concentration Table")],
                [sg.Button("Plot Kinetic Data")],
                [sg.Button("Plot Kinetic Data for Selected Sheets and Species")],
                [sg.Button("Plot Concentration Profiles")],

                [sg.Frame("Settings", layout=[
                    [sg.Text("Concentration unit:"),
                     sg.InputText(key="y_unit", size=(5, 1),default_text=str(data_plotting_settings["y_unit"])),
                     sg.Text("Time unit:"),
                     sg.InputText(key="t_unit", size=(5, 1),default_text=str(data_plotting_settings["t_unit"]))],
                    [sg.Text("Figure size scaler:"),
                     sg.InputText(key="size_scaler", size=(5, 1),default_text=float(data_plotting_settings["size_scaler"])),
                     sg.Text("Max y:"),
                     sg.InputText(key="ylim", size=(5, 1),default_text=str(data_plotting_settings["ylim"]))],
                    [sg.Text("Datapoint size scaler:"),
                     sg.InputText(key="DP_scaler", size=(5, 1),default_text=float(data_plotting_settings["DP_scaler"]))],
                    [sg.Text("Significant figures for initial conc. table:"),
                     sg.InputText(key="SF", size=(5, 1),default_text=float(data_plotting_settings["significant_figs"]))],
                    [sg.Text("Mode for visualising kinetic data:"),sg.Combo(['Scroll', 'Together', 'Separate'], default_value=data_plotting_settings["mode"],
                            key="mode", enable_events=True, readonly=True)],
                ])],

                [sg.Button("Back")]
            ]

            window_inspect_kinetic_data = sg.Window("Inspect Kinetic Data", layout_inspect_kinetic_data)

            while True:
                event_inspect_kinetic_data, values_inspect_kinetic_data = window_inspect_kinetic_data.read()
                if event_inspect_kinetic_data == sg.WIN_CLOSED or event_inspect_kinetic_data == "Back":
                    break
                
                if event_inspect_kinetic_data in ["Generate Initial Concentration Table","Plot Kinetic Data",
                        "Plot Kinetic Data for Selected Sheets and Species","Plot Concentration Profiles"]:
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
                        DP_scaler = float(values_inspect_kinetic_data["DP_scaler"])
                        mode = values_inspect_kinetic_data["mode"]
                        if len(t_unit)==0 or t_unit=='None':
                            t_unit=None
                        size_scaler = float(values_inspect_kinetic_data["size_scaler"])
                        data_plotting_settings["size_scaler"]=size_scaler
                        data_plotting_settings["ylim"]=ylim
                        data_plotting_settings["y_unit"]=y_unit
                        data_plotting_settings["t_unit"]=t_unit
                        data_plotting_settings["mode"]=mode
                        data_plotting_settings["DP_scaler"]=DP_scaler
                        # Check size scaler value
                        if not (0.1 < size_scaler <= 3):
                            raise ValueError("Invalid figure size scaler. Must be between 0.1 and 3.")
                        if not (0 <= DP_scaler <= 5):
                            raise ValueError("Invalid datapoint size scaler. Must be between 0 and 5.")
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
                        data_plotting_settings["significant_figs"]=int(significant_figures)
                    except ValueError as e:
                            sg.popup_error(f"Invalid input: {e}")
                            continue
                if event_inspect_kinetic_data == "Generate Initial Concentration Table":
                    generate_initial_concentration_table(data,size_scaler,y_unit,significant_figures=data_plotting_settings["significant_figs"])
                    
                elif event_inspect_kinetic_data == "Plot Kinetic Data":
                    # Ask the user whether to visualize mass balance using a custom popup
                    layout_visualize_mass_balance = [
                        [sg.Text("Do you want to visualize mass balance for the kinetic data?")],
                        [sg.Button("Yes"), sg.Button("No")]
                    ]

                    window_visualize_mass_balance = sg.Window("Visualize Mass Balance", layout_visualize_mass_balance)

                    event_visualize_mass_balance, values_visualize_mass_balance = window_visualize_mass_balance.read()

                    window_visualize_mass_balance.close()

                    if event_visualize_mass_balance == "Yes":
                        # Use Listbox to create a multi-selection list
                        selected_species = sg.Listbox(
                            values=columns,
                            select_mode=sg.LISTBOX_SELECT_MODE_MULTIPLE,
                            size=(20, min(12, len(fixed_order_species))),
                            key="selected_species_listbox"
                        )

                        layout_select_species = [
                            [sg.Text("Select reaction species for mass balance calculations:")],
                            [selected_species],
                            [sg.Button("OK"), sg.Button("Cancel")]
                        ]

                        window_select_species = sg.Window("Select Reaction Species", layout_select_species)

                        while True:
                            event_select_species, values_select_species = window_select_species.read()

                            if event_select_species == sg.WIN_CLOSED or event_select_species == "Cancel":
                                break
                            elif event_select_species == "OK":
                                selected_species_values = values_select_species["selected_species_listbox"]
                                if not selected_species_values:
                                    sg.popup_error("Please select one or more reaction species.")
                                else:
                                    # Query for stoichiometry numbers
                                    layout_stoichiometry_values = [
                                        [sg.Text(f"Enter stoichiometry values for selected reaction species:")],
                                    ]

                                    for species in selected_species_values:
                                        layout_stoichiometry_values.append(
                                            [sg.Text(f"{species}"), sg.InputText(key=f"stoichiometry_value_{species}", size=(5, 1))])

                                    layout_stoichiometry_values.append([sg.Button("OK"), sg.Button("Cancel")])

                                    window_stoichiometry_values = sg.Window("Enter Stoichiometry Values", layout_stoichiometry_values)

                                    while True:
                                        event_stoichiometry_values, values_stoichiometry_values = window_stoichiometry_values.read()

                                        if event_stoichiometry_values == sg.WIN_CLOSED or event_stoichiometry_values == "Cancel":
                                            break
                                        elif event_stoichiometry_values == "OK":
                                            stoichiometry_values = {}
                                            for species in selected_species_values:
                                                stoichiometry_value = values_stoichiometry_values.get(f"stoichiometry_value_{species}", "")
                                                stoichiometry_values[species] = stoichiometry_value

                                            # Call your plot_kinetic_data function with selected_species_values and stoichiometry_values
                                            plot_kinetic_data(data, stoichiometry_values,ylim=ylim, y_unit=y_unit, t_unit=t_unit, fig_size_scaler=size_scaler,plot_mode=mode,DP_scaler=DP_scaler)
                                            break

                                    window_stoichiometry_values.close()
                                break

                        window_select_species.close()
                    else:
                        # If the user chooses not to visualize mass balance, return None for mass balance and stoichiometry list
                        plot_kinetic_data(data,None,ylim=ylim, y_unit=y_unit, t_unit=t_unit, fig_size_scaler=size_scaler,plot_mode=mode,DP_scaler=DP_scaler)
                elif event_inspect_kinetic_data == "Plot Kinetic Data for Selected Sheets and Species":  # Add this block
                    selected_sheets = values["sheets"]
                    selected_species = values["columns"]
                    if not selected_sheets and not selected_species:
                        sg.popup_error("No sheets nor reaction species have been selected.\n Kinetic data for every reaction species in every sheet will be plotted.")
                    else:
                        if not selected_species:
                            sg.popup_error("No reaction species have been selected.\n Kinetic data for every reaction species in the selected sheets will be plotted.")
                        elif not selected_sheets:
                            sg.popup_error("No sheets have been selected.\n Kinetic data for the selected reaction species in every sheet will be plotted.")
                        # Create a modified kinetic data variable. First, add the dataframes for each selected sheet.
                        data_2=copy.deepcopy(data)
                        data_modified={}
                        for exp in selected_sheets:
                            data_modified[exp]=data_2[exp]
                        # Remove the concentration profiles not included in "selected_species" from the kinetic data.
                        for RM in data_modified.keys():
                            for i, NS in enumerate(data_modified[RM]):
                                if NS not in selected_species and i!=0:
                                    del data_modified[RM][NS]
                        # Ask the user whether to visualize mass balance using a custom popup
                        layout_visualize_mass_balance = [
                            [sg.Text("Do you want to visualize mass balance for the kinetic data?")],
                            [sg.Button("Yes"), sg.Button("No")]
                        ]
                        window_visualize_mass_balance = sg.Window("Visualize Mass Balance", layout_visualize_mass_balance)

                        event_visualize_mass_balance, values_visualize_mass_balance = window_visualize_mass_balance.read()

                        window_visualize_mass_balance.close()

                        if event_visualize_mass_balance == "Yes":
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

                            while True:
                                event_select_species, values_select_species = window_select_species.read()

                                if event_select_species == sg.WIN_CLOSED or event_select_species == "Cancel":
                                    break
                                elif event_select_species == "OK":
                                    selected_species_values = values_select_species["selected_species_listbox"]
                                    if not selected_species_values:
                                        sg.popup_error("Please select one or more reaction species.")
                                    else:
                                        # Query for stoichiometry numbers
                                        layout_stoichiometry_values = [
                                            [sg.Text(f"Enter stoichiometry values for selected reaction species:")],
                                        ]

                                        for species in selected_species_values:
                                            layout_stoichiometry_values.append(
                                                [sg.Text(f"{species}"), sg.InputText(key=f"stoichiometry_value_{species}", size=(5, 1))])

                                        layout_stoichiometry_values.append([sg.Button("OK"), sg.Button("Cancel")])

                                        window_stoichiometry_values = sg.Window("Enter Stoichiometry Values", layout_stoichiometry_values)

                                        while True:
                                            event_stoichiometry_values, values_stoichiometry_values = window_stoichiometry_values.read()

                                            if event_stoichiometry_values == sg.WIN_CLOSED or event_stoichiometry_values == "Cancel":
                                                break
                                            elif event_stoichiometry_values == "OK":
                                                stoichiometry_values = {}
                                                for species in selected_species_values:
                                                    stoichiometry_value = values_stoichiometry_values.get(f"stoichiometry_value_{species}", "")
                                                    stoichiometry_values[species] = stoichiometry_value

                                                # Call your plot_kinetic_data function with selected_species_values and stoichiometry_values
                                                plot_kinetic_data(data_modified, stoichiometry_values,ylim=ylim, y_unit=y_unit, t_unit=t_unit, fig_size_scaler=size_scaler,plot_mode=mode,DP_scaler=DP_scaler)
                                                break

                                        window_stoichiometry_values.close()
                                    break                
                            window_select_species.close()
                        else:
                            # If the user chooses not to visualize mass balance, return None for mass balance and stoichiometry list
                            plot_kinetic_data(data_modified,None,ylim=ylim, y_unit=y_unit, t_unit=t_unit, fig_size_scaler=size_scaler,plot_mode=mode,DP_scaler=DP_scaler)

                elif event_inspect_kinetic_data == "Plot Concentration Profiles":
                    # Create a list of reaction species (assuming fixed_order_species is available)
                    reaction_species_list = fixed_order_species

                    # Use Combo to create a dropdown menu
                    selected_species = sg.Combo(
                        values=reaction_species_list,
                        default_value=reaction_species_list[0] if reaction_species_list else "",
                        size=(20, 1),
                        key="selected_species_dropdown",
                    )

                    layout_select_species = [
                        [sg.Text("Select one reaction species for concentration profiles:")],
                        [selected_species],
                        [sg.Checkbox("Only include selected sheets",key='sheets',default=Plot_concentration_profiles_selected_sheets,enable_events=True)],
                        [sg.Button("OK"), sg.Button("Cancel")]
                    ]

                    window_select_species = sg.Window("Select Reaction Species", layout_select_species)

                    while True:
                        event_select_species, values_select_species = window_select_species.read()
                        if event_select_species == sg.WIN_CLOSED or event_select_species == "Cancel":
                            break
                        elif event_select_species == "OK":
                            Plot_concentration_profiles_selected_sheets=values_select_species['sheets']
                            selected_sheets = values["sheets"]
                            if Plot_concentration_profiles_selected_sheets and len(selected_sheets)>0:
                                data_2=copy.deepcopy(data)
                                data_modified={}
                                for exp in selected_sheets:
                                    data_modified[exp]=data_2[exp]
                            selected_species_value = values_select_species["selected_species_dropdown"]
                            if selected_species_value:
                                if Plot_concentration_profiles_selected_sheets and len(selected_sheets)>0:
                                    plot_concentration_profiles(data_modified, selected_species_value,ylim=ylim, y_unit=y_unit, t_unit=t_unit, fig_size_scaler=size_scaler,DP_scaler=DP_scaler)
                                else:
                                    plot_concentration_profiles(data, selected_species_value,ylim=ylim, y_unit=y_unit, t_unit=t_unit, fig_size_scaler=size_scaler,DP_scaler=DP_scaler)
                            break

                    window_select_species.close()

            window_inspect_kinetic_data.close()
    window.close()

if __name__ == "__main__":
    main()