# Read the SAMtools stat report
with open(input_file, 'r') as file:
    report_text = file.read()

# Extract the relevant statistics using regular expressions
metrics = {
    'Total Sequences': r'raw total sequences:\s+(\d+)',
    'Filtered Sequences': r'filtered sequences:\s+(\d+)',
    'Sequences Mapped': r'reads mapped:\s+(\d+)',
    'Sequences Unmapped': r'reads unmapped:\s+(\d+)',
    'Sequences Duplicated': r'reads duplicated:\s+(\d+)',
    'Reads MQ0': r'reads MQ0:\s+(\d+)',
    'Bases Mapped': r'bases mapped:\s+(\d+)',
    'Bases Trimmed': r'bases trimmed:\s+(\d+)',
    'Mismatches': r'mismatches:\s+(\d+)',
    'Error Rate': r'error rate:\s+([\d.]+)',
    'Average Length': r'average length:\s+(\d+)',
    'Average Quality': r'average quality:\s+([\d.]+)'
}

for metric, pattern in metrics.items():
    metric_data = {}
    match = re.search(pattern, report_text)
    if match:
        value = match.group(1)
        metric_data[metric] = float(value)

    # Plot the metric data
    metric_label = list(metric_data.keys())
    metric_value = list(metric_data.values())

    plt.bar(metric_label, metric_value)
    plt.xlabel('Metrics')
    plt.ylabel('Value')
    plt.title(f'{metric} - SAMtools Statistic')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{metric.lower().replace(" ", "_")}_plot.png'))
    plt.clf()  # Clear the current figure for the next plot

print(f"Plots saved in the '{output_dir}' directory.")
