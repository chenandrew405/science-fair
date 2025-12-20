import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

def create_scatter_plots(df):
    """Create scatter plots for Mean Score vs Matched SNPs."""
    fig, ax = plt.subplots(figsize=(12, 8))

    # Color map for different diseases/drugs
    colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#FFA07A', '#98D8C8', '#F7DC6F']
    markers = ['o', 's', '^', 'D', 'v', 'p']

    for idx, disease in enumerate(df['Disease or Drug']):
        ax.scatter(df.iloc[idx]['Matched SNPs'],
                  df.iloc[idx]['Mean Score'],
                  c=colors[idx],
                  marker=markers[idx],
                  s=300,
                  alpha=0.7,
                  edgecolors='black',
                  linewidth=2,
                  label=disease)

    ax.set_xlabel('Matched SNPs', fontsize=14, fontweight='bold')
    ax.set_ylabel('Mean Score', fontsize=14, fontweight='bold')
    ax.set_title('Mean Risk Score vs Matched SNPs', fontsize=18, fontweight='bold', pad=20)
    ax.legend(loc='upper right', fontsize=10, framealpha=0.9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('data_visualization/scatter_plot.png', dpi=300, bbox_inches='tight')
    print("Scatter plot saved to: data_visualization/scatter_plot.png")
    plt.close()

def create_3d_plots(df):
    """Create 3D scatter plots."""
    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111, projection='3d')

    colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#FFA07A', '#98D8C8', '#F7DC6F']
    markers = ['o', 's', '^', 'D', 'v', 'p']

    for idx, disease in enumerate(df['Disease or Drug']):
        ax.scatter(df.iloc[idx]['Mean Score'],
                  df.iloc[idx]['Matched SNPs'],
                  df.iloc[idx]['Min-Max Normalization'],
                  c=colors[idx],
                  marker=markers[idx],
                  s=400,
                  alpha=0.7,
                  edgecolors='black',
                  linewidth=2,
                  label=disease)

    ax.set_xlabel('Mean Score', fontsize=12, fontweight='bold', labelpad=10)
    ax.set_ylabel('Matched SNPs', fontsize=12, fontweight='bold', labelpad=10)
    ax.set_zlabel('Min-Max Normalization', fontsize=12, fontweight='bold', labelpad=10)
    ax.set_title('3D Visualization of Risk Metrics', fontsize=18, fontweight='bold', pad=20)
    ax.legend(loc='upper left', fontsize=12, framealpha=0.9)

    plt.tight_layout()
    plt.savefig('data_visualization/3d_plot.png', dpi=300, bbox_inches='tight')
    print("3D plot saved to: data_visualization/3d_plot.png")
    plt.close()

def create_bar_graphs(df):
    """Create bar graphs showing Scores vs SNP Ratio."""
    fig, ax = plt.subplots(figsize=(14, 8))

    colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#FFA07A', '#98D8C8', '#F7DC6F']

    diseases = df['Disease or Drug'].tolist()
    mean_scores = df['Mean Score'].tolist()
    matched_snps = df['Matched SNPs'].tolist()

    x = np.arange(len(diseases))
    width = 0.35

    bars1 = ax.bar(x - width/2, mean_scores, width, label='Mean Score',
                   color=colors, alpha=0.8, edgecolor='black', linewidth=1.5)
    bars2 = ax.bar(x + width/2, matched_snps, width, label='Matched SNPs',
                   color=colors, alpha=0.5, edgecolor='black', linewidth=1.5, hatch='//')

    ax.set_xlabel('Disease or Drug', fontsize=14, fontweight='bold')
    ax.set_ylabel('Value', fontsize=14, fontweight='bold')
    ax.set_title('Mean Score vs Matched SNPs by Disease/Drug', fontsize=18, fontweight='bold', pad=20)
    ax.set_xticks(x)
    ax.set_xticklabels(diseases, rotation=45, ha='right', fontsize=11)
    ax.legend(fontsize=12, loc='upper left', framealpha=0.9)
    ax.grid(True, alpha=0.3, axis='y')

    # Add value labels on bars
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax.text(bar.get_x() + bar.get_width()/2., height,
                       f'{height:.2f}',
                       ha='center', va='bottom', fontsize=9, fontweight='bold')

    plt.tight_layout()
    plt.savefig('data_visualization/bar_graph.png', dpi=300, bbox_inches='tight')
    print("Bar graph saved to: data_visualization/bar_graph.png")
    plt.close()

def create_line_graphs(df):
    """Create line graphs showing trends across normalizations."""
    fig, ax = plt.subplots(figsize=(14, 8))

    colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#FFA07A', '#98D8C8', '#F7DC6F']
    markers = ['o', 's', '^', 'D', 'v', 'p']

    norm_cols = ['Mean Score', 'Quantile Normalization', 'Log Normalization', 'Min-Max Normalization']
    x_positions = np.arange(len(norm_cols))

    for idx, disease in enumerate(df['Disease or Drug']):
        values = [df.iloc[idx][col] for col in norm_cols]
        ax.plot(x_positions, values,
               marker=markers[idx],
               color=colors[idx],
               linewidth=2.5,
               markersize=10,
               alpha=0.8,
               label=disease)

    ax.set_xticks(x_positions)
    ax.set_xticklabels(norm_cols, rotation=15, ha='right', fontsize=11)
    ax.set_ylabel('Risk Score Value', fontsize=14, fontweight='bold')
    ax.set_title('Risk Score Across Different Normalization Methods', fontsize=18, fontweight='bold', pad=20)
    ax.legend(loc='upper left', fontsize=11, framealpha=0.9, bbox_to_anchor=(1, 1))
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('data_visualization/line_graph.png', dpi=300, bbox_inches='tight')
    print("Line graph saved to: data_visualization/line_graph.png")
    plt.close()

def create_table():
    """Create a styled table with light blue column headers."""
    # Load the data
    df = pd.read_csv('data_visualization/final_analysis.csv')

    # Rename columns for cleaner display
    df = df.rename(columns={
        'Average Risk Score': 'Mean Score',
        'Average Matched Single Nucleotide Polymorphisms': 'Matched SNPs'
    })

    # Create figure and axis with more width for spacing
    fig, ax = plt.subplots(figsize=(28, 8))
    ax.axis('off')

    # Create the table
    table = ax.table(
        cellText=df.values,
        colLabels=df.columns,
        loc='center',
        cellLoc='center'
    )

    # Styling to match requested format with more spacing
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1, 3)  # Increased vertical scale for more room

    # Flexible column widths - wider for text columns, narrower for numeric
    col_widths = [
        0.11,   # Disease or Drug - wider for names
        0.09,   # Mean Score
        0.09,   # Matched SNPs
        0.07,   # Median
        0.07,   # Mode
        0.07,   # Minimum
        0.07,   # Maximum
        0.07,   # Range
        0.10,   # Quantile Normalization
        0.10,   # Log Normalization
        0.10,   # Min-Max Normalization
        0.07    # Skewness
    ]

    for i, width in enumerate(col_widths):
        for row in range(len(df) + 1):
            table[(row, i)].set_width(width)

    # Header color - light blue
    for (row, col), cell in table.get_celld().items():
        if row == 0:
            cell.set_facecolor("#b0e0e6")  # light blue header
            cell.set_text_props(weight='bold')
            cell.set_height(0.08)  # More height for headers

    # Add title
    plt.title("Genomic Risk Analysis - Final Results", fontsize=22, pad=20)

    # Save the table
    output_file = 'data_visualization/final_analysis_table.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Table saved to: {output_file}")

    # Show the table
    plt.show()

if __name__ == '__main__':
    # Load data with renamed columns
    df = pd.read_csv('data_visualization/final_analysis.csv')
    df = df.rename(columns={
        'Average Risk Score': 'Mean Score',
        'Average Matched Single Nucleotide Polymorphisms': 'Matched SNPs'
    })

    print("Creating all visualizations...")
    print("=" * 60)

    # Create all visualizations
    create_scatter_plots(df)
    create_3d_plots(df)
    create_bar_graphs(df)
    create_line_graphs(df)
    create_table()

    print("=" * 60)
    print("All visualizations created successfully!")
