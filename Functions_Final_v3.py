import numpy as np
import pandas as pd
import os
from itertools import product
from collections.abc import Iterable


# calc possible conc
def allowed_output(value, reaction_vol_nl=20000, drop_size_nl=100, verbose=0):
    # droplet size along with stock conc and total reaction volume restrict number of possible conc to make
    # here we calc all possible conc

    if value['Conc_Values']:
        if isinstance(value['Conc_Stock'], Iterable):
            drop_nums = [i * reaction_vol_nl / (drop_size_nl * value['Conc_Stock'][find_stock(value['Conc_Values'], value['Conc_Stock'], i)[0]]) for i in value['Conc_Values']]
            calculated_concs = value['Conc_Values']
        else:
            drop_nums = [i * reaction_vol_nl / (drop_size_nl * value['Conc_Stock']) for i in value['Conc_Values']]
            calculated_concs = value['Conc_Values']

    else:
        drop_nums = list(range(int((value['Conc_Min'] * reaction_vol_nl) / (drop_size_nl * value['Conc_Stock'])),
                               int((value['Conc_Max'] * reaction_vol_nl) / (drop_size_nl * value['Conc_Stock'])) + 1))

        calculated_concs = [drop_num * value['Conc_Stock'] * drop_size_nl / reaction_vol_nl for drop_num in drop_nums]

    if verbose:
        print('drops :', drop_nums)
        print('volumes :', [i * drop_size_nl for i in drop_nums])
        print('possible_concentrations :', calculated_concs)
    else:
        return calculated_concs, [i * drop_size_nl for i in drop_nums]


# calc executable percentage
def percentage_possible(data, threshold=40):
    lists = list(data.values())

    m = [len(i) for i in data.values()]

    total = np.prod(np.array([len(i) for i in data.values()]))
    possible = 0

    for items in product(*lists):
        if sum(items) <= threshold:
            possible += 1
    
    return (possible/total*100), total

# find best stock for each concentration
def find_stock(conc_values, conc_stocks, this_value):
    num = len(conc_stocks)
    avg = len(conc_values) / float(num)
    out = []
    last = 0.0

    while last < len(conc_values):
        out.append(conc_values[int(last):int(last + avg)])
        last += avg

    for i, value in enumerate(out):
        if this_value in value:
            return i, out

# define random combination generator function_v3.0
def random_combination_generator(concentrations_limits, number_of_combination=100, reaction_vol_nl=10000,
                                 max_nl=None, drop_size_nl=100, check_repeat=True, rounded=2, verbose=0, make_csv=False, return_df=False):
    #  drop size safe
    #  water <0 safe
    #  concentrations_limits is a Dict in this format:
    #  {'name of metabolite': {'Conc_Min': #, 'Conc_Max': #, 'Conc_Values': #, 'Conc_Stock': #, 'Alternatives': #}}

    # generating random combinations
    combinations = []
    data_point = 0
    while data_point < number_of_combination:
        input_data = []
        input_vol = []
        # verbosity
        if (data_point % 10000 == 0) and verbose:
            print(data_point)

        # generation of random input
        for key, value in concentrations_limits.items():
            # Manual Concentration Value Generation
            if value['Conc_Values']:
                # With Alternatives
                if value['Alternatives']:
                    num_alternative = len(value['Alternatives'])
                    choice_alternative = np.random.randint(0, num_alternative)
                    choice_list = [0 for i in range(num_alternative)]
                    choice_list[choice_alternative] = 1

                    choice_conc = np.random.choice(value['Conc_Values'])
                    input_data.append(choice_conc)
                    input_data += choice_list
                    if isinstance(value['Conc_Stock'], Iterable):
                        choice_stock, _ = find_stock(value['Conc_Values'], value['Conc_Stock'], choice_conc)
                        input_vol.append(choice_conc/value['Conc_Stock'][choice_stock]*reaction_vol_nl)
                    else:
                        input_vol.append(choice_conc/value['Conc_Stock']*reaction_vol_nl)

                # Without Alternatives
                else:
                    choice_conc = np.random.choice(value['Conc_Values'])
                    input_data.append(choice_conc)
                    if isinstance(value['Conc_Stock'], Iterable):
                        choice_stock, _ = find_stock(value['Conc_Values'], value['Conc_Stock'], choice_conc)
                        input_vol.append(choice_conc/value['Conc_Stock'][choice_stock]*reaction_vol_nl)
                    else:
                        input_vol.append(choice_conc/value['Conc_Stock']*reaction_vol_nl)

            # Auto Concentration Value Generation
            else:
                # With Alternatives
                if value['Alternatives']:
                    num_alternative = len(value['Alternatives'])
                    choice_alternative = np.random.randint(0, num_alternative)
                    choice_list = [0 for i in range(num_alternative)]
                    choice_list[choice_alternative] = 1

                    drop_num = np.random.randint(round(value['Conc_Min'] * (reaction_vol_nl / drop_size_nl) / value['Conc_Stock']),
                                                 round(value['Conc_Max'] * (reaction_vol_nl / drop_size_nl) / value['Conc_Stock']) + 1)

                    recalculated_conc = drop_num * value['Conc_Stock'] * drop_size_nl / reaction_vol_nl
                    input_data.append(recalculated_conc)
                    input_data += choice_list
                    input_vol.append(recalculated_conc/value['Conc_Stock']*reaction_vol_nl)

                # Without Alternatives
                else:
                    drop_num = np.random.randint(round(value['Conc_Min'] * (reaction_vol_nl / drop_size_nl) / value['Conc_Stock']),
                                                 round(value['Conc_Max'] * (reaction_vol_nl / drop_size_nl) / value['Conc_Stock']) + 1)

                    recalculated_conc = drop_num * value['Conc_Stock'] * drop_size_nl / reaction_vol_nl
                    input_data.append(recalculated_conc)
                    input_vol.append(recalculated_conc/value['Conc_Stock']*reaction_vol_nl)
        
        # Checks
        if check_repetitive and max_nl:
            if input_data not in combinations and sum(input_vol)<= max_nl:
                combinations.append(input_data)
                data_point += 1
        elif check_repetitive and not max_nl:
            if input_data not in combinations:
                combinations.append(input_data)
                data_point += 1
        elif not check_repetitive and max_nl:
            if sum(input_vol)<= max_nl:
                combinations.append(input_data)
                data_point += 1
        else:
            combinations.append(input_data)
            data_point += 1

    # make column name:
    columns_name = []
    for key, value in concentrations_limits.items():
        if not value['Alternatives']:
            columns_name.append(key)
        else:
            columns_name.append(key)
            alternative_name = ['{}_{}'.format(key, i) for i in value['Alternatives']]
            columns_name += alternative_name

    # making csv file
    if make_csv:
        data = pd.DataFrame(np.array(combinations), columns=columns_name)
        data.to_csv('Random_Combination_1.csv', index=False)

    # making dataframe
    if return_df:
        data = pd.DataFrame(np.array(combinations), columns=columns_name)
        return data

    return np.array(combinations)

# transform concentration DataFrame to volume (nanolitre) DataFrame
def concentration_to_volume(concentrations, concentrations_limits, reaction_mixture_vol_nl=10000,
                            fixed_parts={'Lysate': 0.33, 'Saline': 0.1}, round_deg=3):
    # concentrations is a Pandas DataFrame in this format:
    #   {'name of metabolite': concentration}
    # concentrations_limits is a Dict in this format:
    # concentrations_limits (min, max, stock)
    # caution: concentration unit and metabolite name in concentrations and concentrations_limits must be the same

    # make a copy of original dataframe to avoid further change than can affect that
    data = concentrations.copy(deep=True)
    data_all = data.copy(deep=True)
    data = data[[i for i in data.columns if '_' not in i]]
    data *= reaction_mixture_vol_nl

    for metabolite_name, value in concentrations_limits.items():
        if isinstance(value['Conc_Stock'], Iterable):
            print()
            data[metabolite_name] = [round(data[metabolite_name][i] / value['Conc_Stock'][find_stock(value['Conc_Values'], value['Conc_Stock'], data_all[metabolite_name][i])[0]], round_deg) for i in range(len(data[metabolite_name]))]
        else:
            data[metabolite_name] /= value['Conc_Stock']

    # add fix parts
    if fixed_parts:
        for key, value in fixed_parts.items():
            data[key] = reaction_mixture_vol_nl * value

    # add water to reach the reaction_mixture_vol_nl
    data['water'] = reaction_mixture_vol_nl - data.sum(axis=1)

    # for low stock concentration that is not possible to make, raise an error
    # stock conc should set in a way that dont raise this error to avoid further debugging
    if not all(data['water'] > 0): raise Exception("Oops, too concentrated combination!")

    # add alternative
    # make columns name list:
    columns_name = []
    Type_dic = {}
    Stock_dic = {}
    for key, value in concentrations_limits.items():
        if value['Alternatives']:
            columns_name.append(key)
            columns_name.append('{}_Type'.format(key))
            Type_dic[key] = []
        else:
            columns_name.append(key)
        if isinstance(value['Conc_Stock'], Iterable):
            columns_name.append('{}_Stock_Type'.format(key))
            Stock_dic[key] = []

    # Alternatives
    for key in Type_dic.keys():
        data_type = data_all[[i for i in data_all.columns if '{}_'.format(key) in i]]
        for i in data_type.values:
            Type_dic[key].append(concentrations_limits[key]['Alternatives'][np.where(i == 1.0)[0][0]])

    Type_list = list(Type_dic.keys())
    for key in Type_list:
        Type_dic['{}_Type'.format(key)] = Type_dic.pop(key)

    # Stock
    for key in Stock_dic.keys():
        Stock_dic[key] = [concentrations_limits[key]['Conc_Stock'][find_stock(concentrations_limits[key]['Conc_Values'], concentrations_limits[key]['Conc_Stock'], i)[0]] for i in data_all[key]]
         
    Stock_list = list(Stock_dic.keys())
    for key in Stock_list:
        Stock_dic['{}_Stock_Type'.format(key)] = Stock_dic.pop(key)

    data_final = pd.concat([data, pd.DataFrame(Type_dic), pd.DataFrame(Stock_dic)], axis=1)
    return data_final[columns_name + list(fixed_parts.keys()) + ['water']]

# find the first uncompleted day
def day_finder(file, file_format='csv'):
        i = 1
        while True:
            if not os.path.isfile('Day_{}/{}_{}.{}'.format(i, file, i, file_format)):
                return i
            i += 1

# result preprocess
def result_preprocess(day, desired_cols, range=20):
    results = pd.read_csv('Day_{}/Results_{}.csv'.format(day, day))

    # m number pipeline
    data_m = results[desired_cols].iloc[:range, :]
    label_m = results[['yield']].iloc[:range, :]

    # reference, control and specials
    data_specials = results[desired_cols].iloc[range:, :]
    label_specials = results[['yield']].iloc[range:, :]

    return data_m, label_m, data_specials, label_specials

# check to avoid repetitive combinations
def check_repetitive(combination, df_main):
    comparison_df = df_main.merge(combination, indicator=True, how='outer')
    if 'both' in comparison_df._merge.unique():
        return False
    else:
        return True

# main bayesian optimization function
def bayesian_optimization(regressors_list,
                          data, label,
                          concentrations_limits,
                          final_order,
                          reaction_vol_nl=20000, max_nl=13200, drop_size_nl=100,
                          exploitation=1, exploration=1, test_size=100, pool_size=100000, verbose=0, day=1,
                          days_range=[20, 20, 20, 20, 20, 20, 20, 20, 20, 20]):
    
    # first fit training data on our models
    for regressor in regressors_list:
        regressor.fit(data.values, label.values)

    # make random test data
    df_1 = random_combination_generator(concentrations_limits, number_of_combination=pool_size, reaction_vol_nl=reaction_vol_nl,
                                        max_nl=max_nl, drop_size_nl=drop_size_nl, make_csv=False, return_df=True)
    desired_cols = list(df_1.columns)

    df_temp = df_1.copy(deep=True)

    # Upper Confidence Bound
    for index, regressor in enumerate(regressors_list):
        df_1['pred_yield_{}'.format(index)] = regressor.predict(df_temp.values)

    df_1['regressors_std'] = df_1[[str(i) for i in df_1.columns if 'pred_yield' in str(i)]].std(axis=1)
    df_1['mean_vote'] = df_1[[str(i) for i in df_1.columns if 'pred_yield' in str(i)]].mean(axis=1)
    df_1['UCB'] = exploitation * df_1['mean_vote'] + exploration * df_1['regressors_std']
    df_1 = df_1.sort_values(['UCB'], ascending=False)

    # check to don`t make repeated combination but it is not likely
    previous = [pd.read_csv('Day_{}/Concentrations_{}.csv'.format(i, i)).iloc[:days_range[i], :] for i in
                range(1, day + 1)]
    df_main = pd.concat(previous)

    chosen_combinations = pd.DataFrame(columns=desired_cols)
    num = 0
    for i in df_1[desired_cols].values:
        temp_combination = pd.DataFrame([i], columns=desired_cols)
        if check_repetitive(temp_combination, df_main):
            num += 1
            chosen_combinations = pd.concat([chosen_combinations, temp_combination]).reset_index(drop=True)
        if num == test_size:
            break

    return chosen_combinations[final_order]

# ECHO functions
# Part 5: make a dataframe as a 384 well plate for each metabolite
def put_volumes_to_384_wells(volumes_array, starting_well='A1', vertical=False, make_csv=False):
    # vertical and horizontal filling
    # volumes array is concentration_to_volume output
    # volumes array format:
    # a dataframe with columns are component, each row vol of that components
    # this function will output 
    # a list consists of one dataframe for each of metabolite that shows appropriate 384 well plate
    # and one separate dataframe that add well name to volume dataframe 
    if len(volumes_array) > 384: raise ValueError

    all_dataframe = {}
    rows_name = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P']

    if not vertical:
        from_well = rows_name.index(starting_well[0]) * 24 + int(starting_well[1:]) - 1
        # make each metabolite`s dataframe and add it to dict
        for metabolite_name in volumes_array.columns:
            # first make an all zero dataframe
            dataframe = pd.DataFrame(0.0, index=rows_name, columns=range(1, 25))
            # add data one by one in each row
            # (0, 0)--------->(0,23)
            # .......................
            # (15,0)--------->(15,23)
            for index, value in enumerate(volumes_array[metabolite_name]):
                index += from_well
                dataframe.iloc[index // 24, index % 24] = value

            all_dataframe[metabolite_name] = dataframe

        # make named volumes dataframe
        named_volumes = volumes_array.copy(deep=True)
        names = ['{}{}'.format(rows_name[index // 24], index % 24 + 1) for index in named_volumes.index]
        named_volumes['well_name'] = names

    if vertical:
        from_well = rows_name.index(starting_well[0]) + (int(starting_well[1:]) - 1) * 16
        # make each metabolite`s dataframe and add it to dict
        for metabolite_name in volumes_array.columns:
            # first make an all zero dataframe
            dataframe = pd.DataFrame(0.0, index=rows_name, columns=range(1, 25))
            # add data one by one in each column
            # (0, 0)---->-----(0,23)
            # ||||||..........||||||
            # (15,0)---->-----(15,23)
            for index, value in enumerate(volumes_array[metabolite_name]):
                index += from_well
                dataframe.iloc[index % 16, index // 16] = value

            all_dataframe[metabolite_name] = dataframe

        # make named volumes dataframe
        named_volumes = volumes_array.copy(deep=True)
        names = ['{}{}'.format(rows_name[(index + from_well) % 16], (index + from_well) // 16 + 1) for index in
                 named_volumes.index]
        named_volumes['well_name'] = names

    # notice that this function output two value
    return all_dataframe, named_volumes


# make source to destination dataframe for ECHO machine
def source_to_destination(named_volumes, desired_order=None, reset_index=True, check_zero=False):
    # named_volume is second output of put_volumes_to_384_wells function

    # make a separate dataframe for each metabolite that appended in a dict
    # further can be aggregated to one csv file by your desired order
    all_sources = {}
    for metabolite_name in named_volumes.drop(columns=['well_name']):
        transfers = {'Source_Plate_Barcode': [], 'Source_Well': [], 'Destination_Plate_Barcode': [],
                     'Destination_Well': [], 'Transfer_Volume': []}
        for index in range(len(named_volumes)):
            if named_volumes.loc[index, metabolite_name] > 0 or check_zero == False:
                transfers['Source_Plate_Barcode'].append('Plate1')
                transfers['Source_Well'].append('{} well'.format(metabolite_name))
                transfers['Destination_Plate_Barcode'].append('destPlate1')
                transfers['Destination_Well'].append(named_volumes.loc[index, 'well_name'])
                transfers['Transfer_Volume'].append(named_volumes.loc[index, metabolite_name])
        transfers = pd.DataFrame(transfers)

        all_sources[metabolite_name] = transfers

    # aggregate all dataframe
    aggregated = pd.concat(all_sources.values())

    if desired_order:
        aggregated = pd.concat([all_sources[i] for i in desired_order])

    if reset_index:
        aggregated = aggregated.reset_index(drop=True)

    return all_sources, aggregated
