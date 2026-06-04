import numpy as np
import pandas as pd

# import the data
data = pd.read_csv(
    r"C:\tbr\BSc-proj-Alsheghri-2026-tbr\laboratory\16_04_2026_co2_measurements.csv",
    sep=';'
)


data['Timestamp'] = pd.to_datetime(data['Timestamp'])

# start and end for the inlet measurment
start = "2026-04-16 13:00:18"
end   = "2026-04-16 13:07:18"

# extract all the values in the interval
inlet = data[(data["Timestamp"] >= start) & (data["Timestamp"] <= end)]

len(inlet['SCD30_CO2'])
stand = np.std(inlet['SCD30_CO2'])

meann = np.mean(inlet['SCD30_CO2'])
dev = stand/meann * 100 


# oulet concentrations
start_out = "2026-04-16 13:15:18"
end_out = "2026-04-16 14:27:18"
outlet = data[(data["Timestamp"] >= start_out) & (data["Timestamp"] <= end_out)]
len(outlet['SCD30_CO2'])
stand = np.std(outlet['SCD30_CO2'])
meann = np.mean(outlet['SCD30_CO2'])
dev = stand/meann*100