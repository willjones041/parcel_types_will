import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns  # Import Seaborn

# Load the CSV file into a DataFrame
df = pd.read_csv('../filename.csv')

# Get the column headers
column_headers = df.columns.tolist()

# Print the column headers for user to choose from
print("Column Headers:")
for i, header in enumerate(column_headers):
    print(f"{i+1}. {header}")

# User input to select columns to plot
while True:
    try:
        selected_columns = input("Enter the column numbers to plot (separated by comma): ")
        selected_columns = [int(x) - 1 for x in selected_columns.split(',')]

        # Check if two selected columns are provided
        if len(selected_columns) == 2 and all(col in range(len(column_headers)) for col in selected_columns):
            break
        else:
            print("Please select exactly two valid column numbers. Please try again.")
    except ValueError:
        print("Invalid input. Please enter valid column numbers.")

# Set the style for the plot using Seaborn
sns.set(style='darkgrid')  # Use Seaborn's style

# Create the plot
plt.figure(figsize=(10, 6))
plt.plot(df[column_headers[selected_columns[0]]], df[column_headers[selected_columns[1]]], 
         marker='o', linestyle='-', color='b', 
         label=f"{column_headers[selected_columns[0]]} vs {column_headers[selected_columns[1]]}")

# Add legend, labels, and title
plt.legend(fontsize=12)
plt.xlabel(column_headers[selected_columns[0]], fontsize=24)
plt.ylabel(column_headers[selected_columns[1]], fontsize=24)
plt.title(f"{column_headers[selected_columns[1]]} vs {column_headers[selected_columns[0]]}", fontsize=24)
#plt.tick_params(axis='both', labelsize=20)  # Set the font size for the axes ticks
# Add grid and adjust layout
plt.grid(True)
plt.tight_layout()

plt.show()
