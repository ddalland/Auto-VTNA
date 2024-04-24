import numpy as np
import time
import pandas as pd
from polyfit import PolynomRegressor, Constraints
from num2words import num2words
from scipy.optimize import curve_fit
import itertools
import copy
import matplotlib.pyplot as plt
import warnings

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

def VTNA_tT(dt_lists,dNS_lists,order_combination):
        """
        Calculates the numerical integral of the product of the concentrations of normalised species over time \
        raised to their respective order values in order_combination. The numerical integrals are calcualted using the trapezoid rule. 
        Args: 
            time_lists, NS_lists: see explore_orders() 
            order_combination (tuple): the reaction order values used to calculate the normalised time axis. 
        """
        # Make an order combination matrix to match the normalised species 3d array.
        order_array=np.array(order_combination)[:, np.newaxis,np.newaxis]
        # Raise the concentration values of each experiment dataset's normalised species. 
        dNS_lists_updated=np.power(dNS_lists, order_array)
        # Multiply together the 2d normalised species arrays and then multiply with the dt values.
        tT_lists=np.apply_along_axis(lambda x: np.prod(x, axis=0), axis=0, arr=dNS_lists_updated)*dt_lists
        # Obtain the trapezoid rule integrals for each time point using the numpy cumsum function.
        tT_summed=np.cumsum(tT_lists,axis=1)
        # Add 0 to the beginning of each transformed time axis array.
        tT_summed = np.apply_along_axis(lambda x: np.insert(x, 0, 0), axis=1, arr=tT_summed)
        return tT_summed

def data_scaling_y(data,selection):
    """
    Creates a new column in the DataFrames in the kinetic data dictionary "data" consisting of output \
    concentration profiles shifted to start at the same concentration value (only concequential if \
    a reactant rather than product is defined as the output species) and scaled to ensure an average \
    output concentration of 1 across every row of every experiment DataFrame. 
    Args:
            data (dict): Kinetic data as a dictionary of Pandas dataframes.
            VTNA_selection (dict): Dictionary containing information about how to carry out automatic VTNA.
    """
    output_species=selection['output_species']
    VTNA_data=copy.deepcopy(data)    
    species_initial={}
    # Loop through every chosen kinetic profile to collect the first output species datapoint from each experiment.
    if 'RM' in selection.keys():
        experiments=list(selection['RM'].keys())
        for i in selection['RM'].keys():
            species_initial[i]=VTNA_data[i][output_species].to_numpy()[0]
    # If no kinetic profiles have been chosen, find the largest concentration of the output species
    # from all the kinetic profiles in the dataset.
    else: 
        experiments=list(VTNA_data.keys())
        for i in VTNA_data.keys():
            species_initial[i]=VTNA_data[i][output_species].to_numpy()[0]
    #Use largest output concentration value to scale all output concentrations in the VTNA dataset
    max_FOV=max([species_initial[i] for i in list(species_initial.keys())])   
    # Save scaled output concentration values as new column for each reaction profile
    output_values=[]
    for i in experiments:
        shift_value=max_FOV-species_initial[i]
        VTNA_data[str(i)]['output_scaled']=VTNA_data[i][output_species].to_numpy()+shift_value
        output_values=output_values+list(VTNA_data[str(i)]['output_scaled'].to_numpy())
    # Find the average output value across all experiments and use this to scale the 'output_scaled' columns.
    scaling_value=np.mean(output_values)
    for i in experiments:     
        VTNA_data[str(i)]['output_scaled']=np.array([float(i) for i in VTNA_data[i]['output_scaled']])/scaling_value
    return VTNA_data

def on_click(event,VTNA_selection,data,y_unit,size_scaler,fit_metric,deg,constraint,fixed_order_species=None,axis_species=None,show_legend=True):
    """
    Creates a VTNA overlay plot for the kinetic dataset "data" according to the VTNA selection dictionary \
    "VTNA_selection" with order values corresponding to the coordinates of the user's click on a graph or contour \
    plot generated with plot_orders_vs_overlay(). 
    Args:
        event: Contains the x and y positions of the user's click.
        VTNA_selection (dict): Dictionary containing information on the automatic VTNA calculation which generated \
        the overlay score versus order graph or contour plot. 
        data (dict): Kinetic data as a dictionary of Pandas dataframes.
        deg (int): Polynomial degree used in the global fitting prodecure to calculate the overlay score. 
        fit_metric (str): Goodness-of-fit metric chosen to define the overlay score value. Either 'R2', 'RMSE' or 'SE'. \
        otherwise variance is used.
        constraint (str): Defines the constraint of the global fit used by explore_orders() to calculate the overlay score. \
        Can be set to 'monotonic' for monotonic fitting (default), None for ordinary polynomial fitting or 'origin' for \
        linear fit through origin. 
        fixed_order_species (list or None): List containing the normalised reaction species for which \
        time axis should be normalised with a fixed value (as defined in the VTNA selection dictionary). \
        This can be done to lower the dimensionality of the overlay score versus order matrix. 
        axis_species (list): The normalised species of each axis of the overlay score versus order plot. \
        For two normalised species, the x-axis species will appear before the y-axis species. 
        show_legend (bol): Can be set to False to not show the legend of the overlay plot generated. 
    """
    VTNA_selection_copy=copy.deepcopy(VTNA_selection)
    output_species=VTNA_selection_copy['output_species']
    normalised_species_list=list(VTNA_selection_copy['normalised_species'].keys())
    # Create a list of the experiments to include in the fitting procedure. 
    experiments=list(data.keys())
    # If selected experiments RM are defined in the selection dictionary, use these experiment names instead.
    if 'RM' in list(VTNA_selection_copy.keys()):
        experiments=list(VTNA_selection_copy['RM'].keys())
    # Define x and y coordinates as the reaction orders and update the VTNA selection accordingly
    orders=[event.xdata,event.ydata]
    # Define the number of decimals to be included in the x-axis label.
    # Generate data for the second plot based on the click coordinates
    # Define x and y coordinates as the reaction orders and update the VTNA selection accordingly
    siphers=[]
    x_label="$\u03A3"
    counter=0
    if type(fixed_order_species)==str:
        fixed_order_species=[fixed_order_species]
    for i in axis_species:
        VTNA_selection_copy['normalised_species'][i]=orders[counter]
        counter=counter+1
    for j in normalised_species_list:
        q=5
        if str(VTNA_selection_copy['normalised_species'][j])[0]=='-':
            q=q+1
        siphers.append(q)
        x_label=x_label+f"[{j}]^{{{str(VTNA_selection_copy['normalised_species'][j])[:q]}}}"
    x_label=x_label+"\u0394t$."
    # Generate scaled VTNA table for orders.
    tables=VTNA_new(data,VTNA_selection_copy)
    # Scale the output concentration profiles.
    tables_scaled=data_scaling_y(tables,VTNA_selection_copy)
    # Create string to add to plot title if some normalised species have fixed order values
    title_string=f'VTNA overlay plot for {output_species}'
    if fixed_order_species!=None:
        if type(fixed_order_species)==str:
            title_string=title_string+f"\nFixed order for {fixed_order_species}: {VTNA_selection_copy['normalised_species'][fixed_order_species]}"
        else:
            for i in fixed_order_species:
                title_string=title_string+f"\nFixed order for {i}: {VTNA_selection_copy['normalised_species'][i]}. "
    fig2, ax2 = plt.subplots(figsize=(5*size_scaler,4*size_scaler))
    # Check if the output concentration profiles start at the same values for each experiment.
    if 'RM' in list(VTNA_selection.keys()):
        initial_output_values=[data[RM][VTNA_selection['output_species']].to_numpy()[0] for RM in list(VTNA_selection['RM'].keys())]
    else: 
        initial_output_values=[data[RM][VTNA_selection['output_species']].to_numpy()[0] for RM in list(data.keys())]
    if all(elem==initial_output_values[0] for elem in initial_output_values):
        shifted=False
    else:
        shifted=True
    if event.button==3 or event.button==2:
        # Check difference between first and last output values to judge whether fit polynomial should 
        # be monotonically increasing or decreasing.
        delta_output=0
        for RM in experiments:
            if (tables_scaled[RM]['output_scaled'][0]-tables_scaled[RM]['output_scaled'].to_numpy()[-1])>0:
                delta_output=delta_output-1
            else: 
                delta_output=delta_output+1
        fit_direction='dec'
        if delta_output>0:
            fit_direction='inc'
        polyestimator = PolynomRegressor(deg=deg, regularization = None, lam = 0)
        monotone_constraint = Constraints(monotonicity=fit_direction)
        # Plot the time-normalised reaction profiles for each experiment dataset in tables
        Data_points_y=[]
        Data_points_t=[]
        times=[]
        for i in experiments:
            times.append(max([float(j) for j in tables_scaled[i]['tT'].to_numpy()]))
        max_time=max(times)
        for i in experiments:
            # Define the x and y values. 
            x_values=[float(j)/max_time for j in tables_scaled[i]['tT'].to_numpy()]
            Data_points_t=Data_points_t+x_values
            ys=[float(j) for j in tables_scaled[i]['output_scaled'].to_numpy()]
            Data_points_y=Data_points_y+ys
            # Plot the datapoints of the concentration profile 
            ax2.plot(x_values,ys,marker='o',label=i)  
        # Set teh y and x axis labels of the plot. 
        ax2.set_ylabel(f'[{output_species}] / {y_unit}')
        ax2.set_xlabel(x_label)
        # Show the legend if show_legend is set to True. 
        if show_legend:
            legend_1=ax2.legend()
            ax2.add_artist(legend_1)
        y_label=f'Scaled [{output_species}]'
        if shifted:
            y_label='Shifted '+y_label
        ax2.set_ylabel(y_label)
        if constraint=='monotonic':
            # Perform the polynomial fit of "deg" order and evaluate fitted function at the relevant "tT" values       
            Data_points_t_2=np.array(Data_points_t).reshape((-1,1))
            polyestimator.fit(Data_points_t_2, Data_points_y, loss = 'l2', constraints={0: monotone_constraint})       
            pred_mon = polyestimator.predict(Data_points_t_2)
            ts_2=np.linspace(0,1,100).reshape((-1,1))
            pred_mon_fit = polyestimator.predict(ts_2)
            fit_plot,=ax2.plot(ts_2,pred_mon_fit,marker='',label=f'Monotonic {deg}. \norder polynomial fit',linestyle='dotted',linewidth=2,color='black')
            # Call the score_fit function to evalue the goodness-of-fit using the fit measuremetn "fit_metric"
            FQ=score_fit(Data_points_y,pred_mon,fit_metric)
        elif constraint=='origin':
            params = curve_fit(origin_line, Data_points_t, Data_points_y)
            a=params[0]
            y_pred=a*Data_points_t
            FQ=score_fit(Data_points_y,y_pred,fit_metric)
            ts_2=np.linspace(0,max(Data_points_t),100)
            fit_plot,=ax2.plot(ts_2,ts_2*a,marker='',label='linear fit\nvia origin',linestyle='dotted',linewidth=2,color='black')
        else:
            if constraint not in [None,'None']:
                warnings.warn("Ordinary polynomial fitting applied due to unknown constraint input", UserWarning)
            polyfit=np.polyfit(Data_points_t,Data_points_y,deg=deg)
            pred_y=np.polyval(polyfit,Data_points_t)
            tT_2=np.linspace(0,1,100)
            pred_y_fit=np.polyval(polyfit,tT_2)
            FQ=score_fit(Data_points_y,pred_y,fit_metric)
            fit_plot,=ax2.plot(tT_2,pred_y_fit,marker='',label=f'{deg}. order\npolynomial fit',linestyle='dotted',linewidth=2,color='black')
        FQ_round=round(FQ,6)
        plt.title(title_string+f'\nGOF ({fit_metric}): {FQ_round}')
        if show_legend:
            ax2.legend(handles=[fit_plot],loc='lower center')
        plt.tight_layout()
        plt.show()
    elif event.button==1:
        # This corresponds to a left-click.
        # Plot the time-normalised reaction profiles for each experiment dataset in tables
        for i in experiments:
            x_values=[float(j) for j in tables_scaled[i]['tT'].to_numpy()]
            ys=[float(j) for j in tables_scaled[i][output_species].to_numpy()]
            ax2.set_ylabel(f'[{output_species}] / {y_unit}')
            ax2.plot(x_values,ys,marker='o',label=i)
        if show_legend:
            ax2.legend()
        # Define the number of decimals to be included in the x-axis label.
        q,w=5,5
        if str(orders[0])[0]=='-':
            q=1+q
        if str(orders[-1])[0]=='-':
            w=1+w
        # Define the axis labels and the title of the overlay plot. 
        ax2.set_xlabel(x_label)
        y_label=f'Concentration in {y_unit}'
        if shifted:
            y_label='Shifted '+y_label
        ax2.set_ylabel(y_label)
        plt.title(title_string)
        plt.tight_layout()
        plt.show()

def score_fit(y, y_fit, fit_metric):
    # NB: this function was written with the help of OpenAI's ChatGTP. 
    # Find the variance of the residuals of the fit
    FQ = np.var(np.array(y_fit) - np.array(y))
    # Update the goodness-of-fit value according to the chosen fit metric.
    if fit_metric == 'R2':
        # Compute R-squared manually
        SS_res = np.sum((np.array(y) - np.array(y_fit))**2)
        SS_tot = np.sum((np.array(y) - np.mean(np.array(y)))**2)
        FQ = 1 - (SS_res / SS_tot)
    elif fit_metric == 'SE':
        # Compute Standard Error of the Estimate
        FQ = FQ**0.5 / len(y)**0.5
    elif fit_metric == 'RMSE':
        # Compute Root Mean Squared Error
        FQ = np.sqrt(np.mean((np.array(y_fit) - np.array(y))**2))
    return FQ

def VTNA_new(data,VTNA_selection):
    """
    First removes experimental data as specified in "VTNA_selection" from "data". Then \
    calculates the normalised time axis using the trapezoid rule for the concentration \
    profiles of the reaction species listed in "normalised_species" with their respective \
    order values.
    Args: 
        data (dict): Kinetic data as a dictionary of Pandas dataframes.
        VTNA_selection (dict): Dictionary containing information about how to carry out automatic VTNA.
    """
    VTNA_table={}
    #Access one of the experimental dataframes
    example_RM=list(data.keys())[0]
    #Define the time axis key as t
    t=data[example_RM].keys()[0] 
    first_output_values={}
    first_output_values2=[]
    output_species=VTNA_selection['output_species']    
    if 'RM' in list(VTNA_selection.keys()):
        #loop through every selected RM dataframe
        for RM in VTNA_selection['RM'].keys():
            #initialise list TN for normalised time axis
            TN=[float(data[RM][t][0])]
            table=data[RM]
            #Apply omissions of selected datapoints from RM dataframe
            if VTNA_selection["RM"][RM]['omissions']!=None:
                if VTNA_selection["RM"][RM]['omissions'][0]=='range':
                    table=table.drop([c for c in range(len(table)) if c/len(table)*100>VTNA_selection["RM"][RM]['omissions'][1] and c/len(table)*100<VTNA_selection["RM"][RM]['omissions'][2]])
                else: 
                    for j in VTNA_selection["RM"][RM]['omissions']:
                        table=table.drop([int(j)])
            #Update the numbers of each row in case a row has been omitted so that next step works
            table=table.set_index([pd.Index([i for i in range(len(table[t]))])])
            #run through the list of data points for an experiment
            for i in range(len(table[t])-1):
                C=1
                #Define dt as time difference between time point i+1 and the previous one i
                dt=float(table[t][i+1])-float(table[t][i])
                #Apply numerical integration of each normalised species selected in VTNA_selection
                for NS in VTNA_selection["normalised_species"].keys():
                    order=float(VTNA_selection["normalised_species"][NS])
                    C=C*((float(table[NS][i+1])+float(table[NS][i]))/2)**order
                #Add the C trapezoid to the overall normalised time axis
                TN.append(TN[i]+dt*C)
            #Define the running time integral TN as the transformed time axis tT
            table['tT']=TN
            VTNA_table[RM]=table
            # Append the first output concentration value to the first output list and dictionary           
            first_output=float(table[output_species].to_numpy()[0])
            first_output_values[RM]=first_output
            first_output_values2.append(first_output)   
    #Code applied if no RM selections are specified:
    if 'RM' not in list(VTNA_selection.keys()):
        #Update the numbers of each row in case a row has been omitted so that next step works
        for d in data.keys(): 
            TN=[0]
            table=data[d]
            table=table.set_index([pd.Index([i for i in range(len(table[t]))])])
            for i in range(len(table[t])-1):
                #running through the list of data points for an experiment
                C=1
                dt=float(table[t][i+1])-float(table[t][i])
                #Apply numerical integration of the normalised species in VTNA_selection dict
                for NS in VTNA_selection["normalised_species"].keys():
                    order=float(VTNA_selection["normalised_species"][NS])
                    C=C*((float(table[NS][i+1])+float(table[NS][i]))/2)**order
                #Add the dC trapezoid to the overall integrated concentration integral
                TN.append(TN[i]+dt*C)
            table['tT']=TN
            VTNA_table[d]=table
            # Append the first output concentration value to the first output list and dictionary
            first_output=float(table[output_species].to_numpy()[0])
            first_output_values[d]=first_output
            first_output_values2.append(first_output)
    # Find the reaction profile with the largest initial concentration of the output species
    max_FOV=max(first_output_values2)
    # Shift the output reaction profiles to ensure identical starting points.
    for RM in VTNA_table.keys():
        shift_value=max_FOV-first_output_values[RM]
        VTNA_table[RM][output_species]=VTNA_table[RM][output_species].to_numpy()+shift_value
    return VTNA_table

def origin_line(x,a):
    # Function used for performing linear fit through origin
    return x*a

def simple_line(x,a,b):
    # Function used for performing linear fit through origin
    return x*a+b

def VTNA_omissions(data,VTNA_selection,range_mode="Percentage"):
    """
    Uses the 'RM' section of VTNA_selection dictionary to remove experiment DataFrames\
    or selected concentration values within  experiment DataFrames. 
    Args:
            data (dict): Kinetic data as a dictionary of Pandas dataframes.
            VTNA_selection (dict): Dictionary containing information about how to carry out automatic VTNA.
    """
    VTNA_table=copy.deepcopy(data)
    if 'RM' in list(VTNA_selection.keys()):
        VTNA_table={}
        #Access one of the experimental dataframes
        example_RM = next(iter(data.keys()))
        #Define the time axis key as t
        t=data[example_RM].keys()[0] 
        #loop through every selected RM dataframe
        for RM in VTNA_selection['RM'].keys():
            #initialise list TN for normalised time axis
            TN=[float(data[RM][t][0])]
            table=data[RM]
            max_t=table[t].to_numpy()[-1]
            #Apply omissions of selected datapoints from RM dataframe
            if VTNA_selection["RM"][RM]['omissions']!=None:
                if VTNA_selection["RM"][RM]['omissions'][0]=='range':
                    if range_mode=="Percentage":
                        table=table.drop([c for c in range(len(table)) if table[t].to_numpy()[c]
                        /max_t*100>VTNA_selection["RM"][RM]['omissions'][1] and table[t].to_numpy()[c]/max_t*100<=VTNA_selection["RM"][RM]['omissions'][2]])
                    else: 
                        table=table.drop([c for c in range(len(table)) if table[t].to_numpy()[c]
                        >VTNA_selection["RM"][RM]['omissions'][1] and table[t].to_numpy()[c]<=VTNA_selection["RM"][RM]['omissions'][2]])
                else: 
                    for j in VTNA_selection["RM"][RM]['omissions']:
                        table=table.drop([int(j)])
            table=table.set_index([pd.Index([i for i in range(len(table[t]))])])
            VTNA_table[RM]=table
    return VTNA_table