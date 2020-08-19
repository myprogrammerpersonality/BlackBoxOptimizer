import numpy as np
import pandas as pd

# calc possible conc
def allowed_output(conc_limit, reaction_vol_nl = 20000, drop_size_nl = 100, verbose=0):
    # droplet size of ECHO Machines along with stock conc restrict number of possible conc to make
    # here we calc possible conc
    
    drop_nums = list(range(int((conc_limit[0]*reaction_vol_nl)/(drop_size_nl*conc_limit[2])),
                           int((conc_limit[1]*reaction_vol_nl)/(drop_size_nl*conc_limit[2]))+1))

    calculated_concs = [drop_num * conc_limit[2] * drop_size_nl / reaction_vol_nl for drop_num in drop_nums]
    if verbose:
        print('drops :', drop_nums)
        print('volumes :', [i * drop_size_nl for i in drop_nums])
        print('possible_concentrations :', calculated_concs)
    else:
        return drop_nums, calculated_concs


# Part 2: define random combination generator funtion
def random_combination_generator(concentrations_limits, number_of_combination=100, reaction_vol_nl=10000,
                                 check_max=True, max_nl=10000, drop_size_nl=25, rounded=2, verbose=0, make_csv=False, return_df=False):
    #  drop size safe
    #  water <0 safe
    #  concentrations_limits is a Dict in this format:
    #  {'name of metabolite': (min, max, stock)}
    #  make this list for checking max vol
    
    stocks_vol = [reaction_vol_nl/i[2] for i in concentrations_limits.values()]
    
    combinations = []
    data_point = 0
    while data_point < number_of_combination:
        input_data = []
        # verbosity
        if (data_point % 10000 == 0) and verbose:
            print(data_point)
        # generation of random input
        for key, value in concentrations_limits.items():
            # these two line make output concentartions, safe for drop_size nl droplet size of ECHO Machines
            drop_num = np.random.randint(round(value[0]*(reaction_vol_nl/drop_size_nl)/value[2]),
                                         round(value[1]*(reaction_vol_nl/drop_size_nl)/value[2])+1)

            recalculated_conc = drop_num * value[2] * drop_size_nl / reaction_vol_nl        
            input_data.append(recalculated_conc)
        
        if check_max:
            if sum([i*j for i, j in zip(input_data, stocks_vol)]) <= max_nl:
                # appending to other combination
                combinations.append(input_data)
                data_point += 1
            else:
                pass
        else:
            # appending to other combination
            combinations.append(input_data)
            data_point += 1
    
    # making csv file
    if make_csv:
        data = pd.DataFrame(np.array(combinations), columns=concentrations_limits.keys())
        data.to_csv('Random_Combination_1.csv',index=False)
        
    # making dataframe
    if return_df:
        data = pd.DataFrame(np.array(combinations), columns=concentrations_limits.keys())
        return data

    return np.array(combinations)

# transform concentration DataFrame to volume (nanolitre) DataFrame
def concentration_to_volume(concentrations, concentrations_limits, reaction_mixture_vol_nl=10000,
                            fixed_parts = {'Lysate':0.33, 'Saline':0.1}, round_3=True):
    ## concentrations is a Pandas DataFrame in this format:
    #   {'name of metabolite': concentration}
    ## concentrations_limits is a Dict in this format:
    # concentrations_limits (min, max, stock)
    ## caution: concentration unit and metabolite name in concentrations and concentrations_limits must be the same
    
    # make a copy of original dataframe to avoid further change to affect that
    data = concentrations.copy(deep=True)
    
    data *= reaction_mixture_vol_nl
    for metabolite_name, value in  concentrations_limits.items():
        stock_conc = value[2]
        if round_3:
            data[metabolite_name] = round(data[metabolite_name] / stock_conc, 3)    
        else:
            data[metabolite_name] /= stock_conc
    #
    ### data['total_except_water'] = data.sum(axis=1)
    #
    
    # add fix parts
    if fixed_parts:
        for key, value in fixed_parts.items():
            data[key] = reaction_mixture_vol_nl * value
    
    # add water to reach the reaction_mixture_vol_nl
    data['water'] = reaction_mixture_vol_nl - data.sum(axis=1)
    
    # for low stock concentration that is not possible to make, raise an error
    # stock conc should set in a way that dont raise this error to avoid further debugging
    if not all(data['water']>0): raise Exception("Oops, too concentrated combination!")
        
    # return vol in nl
    return data


# Part 5: make a dataframe as a 384 well plate for each metabolite
def put_volumes_to_384_wells(volumes_array, starting_well='A1', vertical=False, make_csv=False):
    # vertical and horizontal filling
    # volumes array is concentration_to_volume output
    # volumes array format:
    # a dataframe with columns are component, each row vol of that components
    # this function will output 
    # a list consists of one dataframe for each of metabolite that shows appropriate 384 well plate
    # and one separate dataframe that add well name to volume dataframe 
    if  len(volumes_array) > 384: raise ValueError
    
    all_dataframe = {}
    rows_name = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P']
    
    if not vertical:
        from_well = rows_name.index(starting_well[0])*24 + int(starting_well[1:]) - 1
        # make each meatabolite`s dataframe and add it to dict
        for metabolite_name in volumes_array.columns:
            # first make an all zero dataframe
            dataframe = pd.DataFrame(0.0, index=rows_name, columns=range(1, 25))
            # add data one by one in each row
            # (0, 0)--------->(0,23)
            # .......................
            # (15,0)--------->(15,23)
            for index, value in enumerate(volumes_array[metabolite_name]):
                index += from_well
                dataframe.iloc[index//24, index%24] = value
        
            all_dataframe[metabolite_name]= dataframe
    
        # make named volumes dataframe
        named_volumes = volumes_array.copy(deep=True)
        names = ['{}{}'.format(rows_name[index//24], index%24+1) for index in named_volumes.index]
        named_volumes['well_name'] = names
    
    if vertical:
        from_well = rows_name.index(starting_well[0]) + (int(starting_well[1:])-1)*16
        # make each meatabolite`s dataframe and add it to dict
        for metabolite_name in volumes_array.columns:
            # first make an all zero dataframe
            dataframe = pd.DataFrame(0.0, index=rows_name, columns=range(1, 25))
            # add data one by one in each column
            # (0, 0)---->-----(0,23)
            # ||||||..........||||||
            # (15,0)---->-----(15,23)
            for index, value in enumerate(volumes_array[metabolite_name]):
                index += from_well
                dataframe.iloc[index%16, index//16] = value
        
            all_dataframe[metabolite_name]= dataframe
    
        # make named volumes dataframe
        named_volumes = volumes_array.copy(deep=True)
        names = ['{}{}'.format(rows_name[(index+from_well)%16], (index+from_well)//16+1) for index in named_volumes.index]
        named_volumes['well_name'] = names
    
    # notice that this function output two value
    return all_dataframe, named_volumes


# make source to destination dataframe for ECHO machine
def source_to_destination(named_volumes, desired_order=None, reset_index=True, check_zero=False):
    # named_volume is second output of put_volumes_to_384_wells function
    
    # make a separate dataframe for each matabolite that appended in a dict
    # further can be aggregated to one csv file by your desired order
    all_sources = {}
    for metabolite_name in named_volumes.drop(columns=['well_name']):
        transfers = {'Source_Plate_Barcode':[], 'Source_Well':[], 'Destination_Plate_Barcode':[], 'Destination_Well':[], 'Transfer_Volume':[]}
        for index in range(len(named_volumes)):
            if named_volumes.loc[index, metabolite_name] > 0 or check_zero == False:
                transfers['Source_Plate_Barcode'].append('Plate1')
                transfers['Source_Well'].append('{} well'.format(metabolite_name))
                transfers['Destination_Plate_Barcode'].append('destPlate1')
                transfers['Destination_Well'].append(named_volumes.loc[index, 'well_name'])
                transfers['Transfer_Volume'].append(named_volumes.loc[index, metabolite_name])
        transfers = pd.DataFrame(transfers)
        
        all_sources[metabolite_name]=transfers
    
    # aggregate all dataframe
    aggregated = pd.concat(all_sources.values())
    
    if desired_order:
        aggregated = pd.concat([all_sources[i] for i in desired_order])
    
    if reset_index:
        aggregated = aggregated.reset_index(drop=True)
    
    return all_sources, aggregated


def bayesian_optimization(regressors_list,
                          data,  label,
                          concentrations_limits,
                          desired_cols, fixed_col_value, final_order,
                          exploitation=1, exploration=1, test_size=100, pool_size=100000, verbose=0):
    
    
    # first fit
    for regressor in regressors_list:
        regressor.fit(data.values, label.values)
        
    # make i th dataset
    df_1 = random_combination_generator(concentrations_limits, number_of_combination=pool_size, reaction_vol_nl=20000,
                                 check_max=True, max_nl=13200, drop_size_nl=100, make_csv=False, return_df=True)
 
    df_1 = df_1[desired_cols]
    
    df_temp = df_1.copy(deep=True)
    
    #upper Confidence Bound
    for index, regressor in enumerate(regressors_list):
        df_1['pred_yield_{}'.format(index)] = regressor.predict(df_temp.values)
        
    df_1['regressors_std'] = df_1[[str(i) for i in df_1.columns if 'pred_yield' in str(i)]].std(axis=1)
    df_1['mean_vote'] = df_1[[str(i) for i in df_1.columns if 'pred_yield' in str(i)]].mean(axis=1)
    df_1['UCB'] = exploitation * df_1['mean_vote']+ exploration * df_1['regressors_std'] 
    df_1 = df_1.sort_values(['UCB'], ascending=False)
    
    chosen_combinations = df_1.iloc[0:test_size,0:len(desired_cols)].values
    chosen_combinations = pd.DataFrame(chosen_combinations, columns=desired_cols)
    
    for key, value in fixed_col_value.items():
        chosen_combinations[key] = value[0]

    chosen_combinations = chosen_combinations[final_order]
    
    return chosen_combinations

# result preprocess

def result_preprocess(day, desired_cols, ranges=(20, 24)):
    
    results = pd.read_csv('Day_{}/Results_{}.csv'.format(day, day))

    # 20 number pipeline
    data_20 = results[desired_cols].iloc[:ranges[0],:-1]
    label_20 = results[desired_cols].iloc[:ranges[0],-1:]

    # ref_excel, jove, ours , control
    data_specials = results[desired_cols].iloc[ranges[0]:ranges[1],:-1]
    label_specials = results[desired_cols].iloc[ranges[0]:ranges[1],-1:]
    
    return data_20, label_20, data_specials, label_specials


# find the first uncomplete day
def day_finder(file, file_format='csv'):
    for i in range(1, 12):
        if not os.path.isfile('Day_{}/{}_{}.{}'.format(i, file, i, file_format)):
            return i
    return 0

def process_limits(concentrations_limits):
    fixed = {}
    changing = []
    for key, value in concentrations_limits.items():
        if value[0] == value[1]:
            fixed[key] = value[0]
        else:
            changing.append(key)

    return fixed, changing, list(concentrations_limits.keys())
