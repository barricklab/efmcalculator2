import json
import numpy as np
from scipy.interpolate import BSpline

# Load model data from JSON
with open('gam_model.json', 'r') as file:
    model_data = json.load(file)

# Extract model parameters from loaded JSON
coefficients = np.array(model_data['coefficients'])
RBP_Length_knots = np.array(model_data['RBP_Length_knots'], dtype=float)  # Ensure floats
TBD_length_knots = np.array(model_data['TBD_length_knots'], dtype=float)  # Ensure floats

# Define prediction function in Python
def predict_gam_model(RBP_Length, TBD_length):
    # Implement prediction logic using coefficients and knots
    linear_pred = np.dot(coefficients, [1, RBP_Length, TBD_length, RBP_Length**2, TBD_length**2])  # Example linear prediction
    spline_pred = BSpline(RBP_Length_knots, np.zeros_like(RBP_Length_knots), 3)(TBD_length)
    return linear_pred + spline_pred

# Example usage
predicted_rate = predict_gam_model(6, 50)
print(f"Predicted Mutation Rate: {predicted_rate}")
