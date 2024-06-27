import pandas as pd
from sklearn.model_selection import train_test_split
from pygam import LinearGAM, s

#load csv file
data_of_interest = pd.read_csv('new.csv')
print(data_of_interest.columns)


#separating df
X = data_of_interest[['RBP_Length', 'TBD_length']]
y = data_of_interest['Mutation_Rate']
weights = data_of_interest['weight']

#gam model fit with cr splines and regularization b/c small data
gam_model = LinearGAM(s(0, n_splines=4, spline_order=3) +  s(1, n_splines=4, spline_order=3)) 
gam_model.fit(X, y, weights=weights
              )

#data not included
new_data = pd.DataFrame({'RBP_Length': [6], 'TBD_length': [50]})


#prediction site
predictions = gam_model.predict(new_data)
predictions_abs = abs(predictions)

print(gam_model.summary())
print("Mutation Rate:", predictions_abs)

