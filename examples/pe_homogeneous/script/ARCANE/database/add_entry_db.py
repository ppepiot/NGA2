import sqlite3 as sl
import pandas as pd

# Connect to the database
con = sl.connect('species_tpf_data.db')

# Add the missing data here
data_dict = {}
data_dict['name'] = ["IC17H34"]
data_dict['liquid_density'] = [793.0]
data_dict['boiling_temperature'] = [513.2]
data_dict['critical_temperature'] = [692.0]
data_dict['critical_pressure'] = [15.7]
data_dict['critical_volume'] = [890e-6]
data_dict['viscosity_params_3rd_order'] = [-8.8790]
data_dict['viscosity_params_2nd_order'] = [1.4480E+03]
data_dict['viscosity_params_1rst_order'] = [1.8074E-02]
data_dict['viscosity_params_zero_order'] = [-1.4953E-05]

# Add the new data in the database
new_db_data = pd.DataFrame.from_dict(data_dict)
new_db_data.to_sql('SPECIES_TPF',con, if_exists='append',index=False)

# Delete duplicates
with con:
    con.execute("""
                DELETE FROM SPECIES_TPF
                WHERE rowid NOT IN (
                SELECT MIN(rowid)
                FROM SPECIES_TPF
                GROUP BY name
                )
                """)

#con.execute("DELETE FROM SPECIES_TPF WHERE name == 'IC17H34'")

