import os
import pysam

# Prompt the user for the input BAM file
bam_filename = input("Enter the input BAM file name: ")

# Check if file exists
if not os.path.isfile(bam_filename):
    print(f"File '{bam_filename}' does not exist.")
    exit(1)

bam_file = pysam.AlignmentFile(bam_filename, "rb")

# Check if BAM file is indexed
try:
    bam_file.check_index()
    print(f"BAM file '{bam_filename}' is indexed.")
except ValueError:
    print(f"BAM file '{bam_filename}' is not indexed.")
    bam_file.close()
    exit(1)

# Check if BAM file is sorted
header_dict = bam_file.header.to_dict()
if "HD" in header_dict and "SO" in header_dict["HD"]:
    if header_dict["HD"]["SO"] == "coordinate":
        print(f"BAM file '{bam_filename}' is sorted.")
    else:
        print(f"BAM file '{bam_filename}' is not sorted.")
        bam_file.close()
        exit(1)
else:
    print("Unable to determine if BAM file is sorted.")
    bam_file.close()
    exit(1)

# Split the filename and extension
input_filename, _ = os.path.splitext(os.path.basename(bam_filename))

# Prompt the user for the minimum alignment length
while True:
    min_length = input("Enter the minimum alignment length: ")
    try:
        min_length = int(min_length)
        break
    except ValueError:
        print("Invalid input. Please enter an integer.")

# Output BAM file
out_filename = f"{input_filename}_filtered_{min_length}.bam"
out_bam = pysam.AlignmentFile(out_filename, "wb", header=bam_file.header)

# Filter reads based on alignment length
for read in bam_file:
    if read.query_alignment_length >= min_length:
        out_bam.write(read)

# Close the BAM files
bam_file.close()
out_bam.close()

print(f"Filtered BAM file saved as {out_filename}.")

