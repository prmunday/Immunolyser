# from app import utils
# from constants import *
# import app

# if __name__ == "__main__":
#     # Replace this with relevant test cases for your method
#     result = utils.generateBindingPredictions("6fa3cf66-72a7-4e9d-9165-a4a46cfd38c0", "HLA-A0201", Class_One_Predictors.MHCflurry)
#     print(result)

from ProtPeptigram.DataProcessor import PeptideDataProcessor
from ProtPeptigram.viz import ImmunoViz

# Initialize data processor
processor = PeptideDataProcessor()

# Load data
processor.load_peaks_data("data/peptides.csv")
processor.load_protein_sequences("data/proteome.fasta")

# Process data
formatted_data = processor.filter_and_format_data(
    filter_contaminants=True,
    intensity_threshold=1000,
    min_samples=2
)

# Create visualizations
viz = ImmunoViz(formatted_data)
fig, _ = viz.plot_peptigram(
    protein_ids=["P20152", "P32261"],
    group_by="Sample",
    color_by="protein",
    title="HLA Peptide Visualization"
)

# Save visualization
fig.savefig("protein_visualization.png", dpi=300, bbox_inches="tight")