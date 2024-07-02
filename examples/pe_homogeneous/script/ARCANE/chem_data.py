"""File for the storage of non-computable chemical properties
from external databases

Author: QC (2020/03/24)
Modification : coupled with actual external database , TL (2021/07/29)
"""
import sqlite3 as sl
import pandas as pd
import os
import ARCANE.display as display

logger = display.Logger()
logger.set_log('logCompute')


def get_tpf_data(species_name):
    """Retrieves the two-phase flow block for species calculation

    Parameters
    ----------
    species_name : str
        name of the species

    Returns
    -------
    species_tpf_dict : dict
        dictionary containing the data

    """

    species_tpf_dict = {}

    # Connect to external database
    file = 'file:'+ os.path.dirname(os.path.abspath(__file__)) + '/database/species_tpf_data.db?model=rw'
    try:
        con = sl.connect(file, uri=True)
    except:
        logger.error("Could not connect to the species database")
        quit()

    # Get the informations stored in the external database
    stored_data = pd.read_sql("SELECT * FROM SPECIES_TPF",con)

    # Fuels that are accepted
    species = stored_data['name'].to_list()

    # Check if species is in the database and get it if it is
    query = 'name=="'+species_name+'"'
    species_data = stored_data.query(query)
    if species_data.empty:
        logger.info('Species ' + species_name + ' is not available.\n' +
                    'If you want to compute your own species, you should write directly into the sources of ARCANE, '
                    'using yaws book (see chemistry website)')
        logger.info('Species you can choose now are :')
        for f in species:
            logger.info(f)
    else:
        species_tpf_dict['liquid_density'] = species_data.iloc[0]['liquid_density']
        species_tpf_dict['boiling_temperature'] = species_data.iloc[0]['boiling_temperature']
        species_tpf_dict['critical_temperature'] = species_data.iloc[0]['critical_temperature']
        species_tpf_dict['critical_pressure'] = species_data.iloc[0]['critical_pressure']
        species_tpf_dict['critical_volume'] = species_data.iloc[0]['critical_volume']
        species_tpf_dict['viscosity_params'] = [species_data.iloc[0]['viscosity_params_3rd_order'],
                                                species_data.iloc[0]['viscosity_params_2nd_order'],
                                                species_data.iloc[0]['viscosity_params_1rst_order'],
                                                species_data.iloc[0]['viscosity_params_zero_order']]

    return species_tpf_dict
