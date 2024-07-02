import sqlite3 as sl
import pandas as pd
con = sl.connect('species_tpf_data.db')

with con :
        con.execute("""
        CREATE TABLE SPECIES_TPF(
        name TEXT,
        liquid_density FLOAT,
        boiling_temperature FLOAT,
        critical_temperature FLOAT,
        critical_pressure FLOAT,
        critical_volume FLOAT,
        viscosity_params_3rd_order FLOAT,
        viscosity_params_2nd_order FLOAT,
        viscosity_params_1rst_order FLOAT,
        viscosity_params_zero_order FLOAT);
        """)

sql = 'INSERT INTO SPECIES_TPF (name,liquid_density ,boiling_temperature,critical_temperature,critical_pressure,critical_volume,viscosity_params_3rd_order,viscosity_params_2nd_order, viscosity_params_1rst_order, viscosity_params_zero_order) values(?,?,?,?,?,?,?,?,?,?)'

data = [ ("O2",1141.0, 90.2, 154.58, 50.43, 90e-6, -2.1354, 8.6824E+01, 1.0541E-02, -6.0681E-05),
         ("H2", 70.99, 20.1, 33.18, 13.0, 33e-6, -2.9289, 1.4778E+01, 2.8833E-02, -6.3557E-04),
         ("CH4", 422.0, 111.7, 190.56, 45.99, 100e-6, -7.3801, 3.1925E+02, 4.7934E-02, -1.4120E-04),
         ("C2H2", 567.4, 169.35, 282.34, 50.41, 130e-6, -0.0709, 2.8381E+01, -4.6617E-03, 3.1151E-06),
         ("C3H8", 493.0, 231.1, 369.8, 42.5, 200e-6, -3.1759, 2.9712E+02, 9.5453E-03, -1.8781E-05),
         ("C4H10", 573.0, 272.7, 425.1, 38.0, 260e-6, -1.8077, 2.5893E+02, 3.0205E-03, -8.6441E-06),
         ("C7H8", 862.1, 383.8, 593.0, 41.0, 316.00e-6, -1.9879, 5.0806E+02, 1.2152E-03, -2.7318E-06),
         ("MCYC6", 770.0, 374.09, 572.19, 34.71, 368.00e-6, -1.9879, 5.0806E+02, 1.2152E-03, -2.7318E-06),
         ("NC7H16", 684.0, 371.57, 540.20, 27.40, 428.0e-6, -5.2620, 7.4043E+02, 1.2279E-02, -1.4466E-05),
         ("XYLENE", 777.0, 417.59, 630.30, 37.32, 370e-6, -9.9569, 1.5752E+03, 2.0516E-02, -1.7208E-05),
         ("IC8H18", 690.0, 372.2, 543.9, 25.7, 469.7e-6, -15.0420, 2.0021E+03, 3.71E-02, -3.4486E-05),
         ("DECALIN", 896.0, 463.2, 645.0, 20.80, 480, -3.2074, 8.4444E+02, 3.7053E-03, -3.8865E-06),
         ("NC10H22", 730.0, 447.3, 617.7, 20.03, 620e-6, -6.0716, 1.0177E+03, 1.2247E-02, -1.1892E-05),
         ("C11H10", 1000.0, 514.70, 772.0, 36.0, 470e-6, -8.1601, 1.5562E+03, 1.4221E-02, -1.0281E-05),
         ('NC12H26', 750.0, 489.5, 658.0, 18.0, 750e-6, -7.0687, 1.2530E+03, 1.3735E-02, -1.2215E-05),
         ("IC12H26", 750.0, 451.0, 651.57, 18.51, 680e-6, -7.4683, 1.1508E+03, 1.5874E-02, -1.4233E-05),
         ("IC16H34", 793.0, 513.2, 692.0, 15.7, 890e-6, -8.8790, 1.4480E+03, 1.8074E-02, -1.4953E-05),
         ("NC16H34", 793.0, 554.0, 722.0, 14.0, 0.00103, -8.1894, 1557.1, 0.01527, -1.2371E-05),
         ("TMBENZ", 900.0, 470.0, 676.0, 29.0, 487e-6, -5.1984, 9.7976e+2, 8.8146e-3, -7.5668e-6 ),
         ("TETRALIN", 970.0, 480.0, 720.0, 37.0, 408e-6, -6.3710, 1274.0, 1.0494e-2, -8.1163e-6 ),
         ("C10H7CH3", 1020.0, 514.0, 764.0, 32.9, 466e-6, -9.3427, 1668.9, 1.7940e-2, -1.4115e-5 )]

with con :
        con.executemany(sql,data)

with con:
        con.execute("""
                    DELETE FROM SPECIES_TPF
                    WHERE rowid NOT IN (
                    SELECT MIN(rowid)
                    FROM SPECIES_TPF
                    GROUP BY name
                    )
                    """)

        data = con.execute("""
                           SELECT * FROM SPECIES_TPF
                           """)
        for row in data:
                print(row)

df = pd.read_sql("""
                 SELECT * FROM SPECIES_TPF
                 """,con)

df.to_csv('species_tpf_data.csv',index=False)
