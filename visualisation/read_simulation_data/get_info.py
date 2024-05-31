import h5py

def get_info(folder):
    def extract_group(group):
        group_dict = {}
        for key, item in group.items():
            # Extracting the value for each item
            value = item[()]

            # Extracting the first attribute
            if item.attrs:
                for attr_name in item.attrs:
                    unit_attr = item.attrs[attr_name]
                    unit = unit_attr[0].decode('utf-8') if isinstance(unit_attr[0], bytes) else unit_attr[0]
                    break
            else:
                unit = 'unknown'

            group_dict[key] = [value, unit]
        
        return group_dict

    with h5py.File(f'{folder}info.h5', 'r') as file:
        constants_dict = extract_group(file['constants'])
        global_params_dict = extract_group(file['global_parameters'])

    with h5py.File(f'{folder}snap0_0.h5') as file:
        grid_info_dict = extract_group(file['grid_info'])   

    return constants_dict, global_params_dict, grid_info_dict