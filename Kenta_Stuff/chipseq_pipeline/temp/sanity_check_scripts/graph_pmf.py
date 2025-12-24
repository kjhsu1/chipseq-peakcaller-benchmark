"""
Imports
"""
import argparse
import pandas as pd
import matplotlib.pyplot as plt

"""
Argument Parsing
"""
parser = argparse.ArgumentParser(
    description='Graph the PMF of the given pmf+variance CSV'
)

parser.add_argument('--pmf_var_csv', type=str, required=True,
    help='Path to pmf+variance CSV')

args = parser.parse_args()

"""
Main Logic
"""
def split_by_restart(df, key='bin_idx'):
    """Split the DataFrame every time bin_idx restarts (e.g., goes from high to 0)."""
    splits = []
    current = []
    
    prev = -1
    for _, row in df.iterrows():
        if row[key] == 0 and prev != -1:
            # Restart detected
            splits.append(pd.DataFrame(current))
            current = []
        current.append(row)
        prev = row[key]
    
    # Append the last chunk
    if current:
        splits.append(pd.DataFrame(current))
    
    return splits

def plot_pmf_groups(df_groups):
    for i, group in enumerate(df_groups):
        bins = group['bin_idx']
        pmf = group['pmf']
        
        plt.figure(figsize=(10, 5))
        plt.plot(bins, pmf, marker='o', linestyle='-', color='blue')
        plt.title(f"Probability Mass Function (Chromosome {i+1})")
        plt.xlabel("Bin Index")
        plt.ylabel("PMF Value")
        plt.grid(True)
        plt.tight_layout()
        plt.show()

def main():
    df = pd.read_csv(args.pmf_var_csv)
    df_groups = split_by_restart(df)
    plot_pmf_groups(df_groups)

# Run the program
main()
