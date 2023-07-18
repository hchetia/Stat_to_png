#VisEase-v1.0
#Author- Hasnahana Chetia
#Date 17-07-2023

import os
import matplotlib.pyplot as plt
import pysam

# Load the BAM file
bam_filename = "sorted_7400_1740_filtered_250.sorted.bam"
bam_file = pysam.AlignmentFile(bam_filename, "rb")

# Split the filename and extension
filename, _ = os.path.splitext(bam_filename)

# Compute the alignment lengths
alignment_lengths = []
for read in bam_file:
    if not read.is_unmapped:
        alignment_lengths.append(read.query_alignment_length)

print(f"Collected {len(alignment_lengths)} alignment lengths.")

# Check if there are any alignment lengths
if alignment_lengths:
    # Plot the histogram
    plt.figure(figsize=(10, 6))
    plt.hist(alignment_lengths, bins=100, color='skyblue', edgecolor='black')
    plt.xlabel('Alignment Length')
    plt.ylabel('Frequency')
    plt.title(f'Distribution of Alignment Lengths ({filename})')

    # Save the plot
    output_file = f"{filename}_alignment_lengths.png"
    plt.savefig(output_file)
    print(f"Plot saved as {output_file}.")
else:
    print("No alignment lengths collected, skipping plot generation.")

