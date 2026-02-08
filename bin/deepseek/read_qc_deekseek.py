import argparse
from openai import OpenAI
import re

def parse_nanostats_file(file_path):
    """
    Parse the NanoStats file and extract key sequencing metrics.
    
    Parameters:
    file_path (str): Path to the NanoStats text file.
    
    Returns:
    dict: A dictionary containing the parsed sequencing metrics.
    """
    metrics = {}

    with open(file_path, 'r') as file:
        content = file.read()

    # Regular expressions to match relevant metrics
    def safe_search(pattern, content, default=None):
        match = re.search(pattern, content)
        if match:
            return float(match.group(1)) if default is None else match.group(1)
        return default

    # Extract the data using safe_search
    metrics['mean_read_length'] = safe_search(r"Mean read length:\s*(\d+\.\d+)", content)
    metrics['mean_read_quality'] = safe_search(r"Mean read quality:\s*(\d+\.\d+)", content)
    metrics['median_read_length'] = safe_search(r"Median read length:\s*(\d+\.\d+)", content)
    metrics['median_read_quality'] = safe_search(r"Median read quality:\s*(\d+\.\d+)", content)
    metrics['number_of_reads'] = safe_search(r"Number of reads:\s*(\d+)", content, default=0)
    metrics['read_length_n50'] = safe_search(r"Read length N50:\s*(\d+\.\d+)", content)
    metrics['stdev_read_length'] = safe_search(r"STDEV read length:\s*(\d+\.\d+)", content)
    metrics['total_bases'] = safe_search(r"Total bases:\s*(\d+\.\d+)", content)

    # Parsing quality cutoff information
    quality_cutoffs = re.findall(r">Q(\d+):\s*(\d+)\s*\((\d+\.\d+)%\)\s*(\d+\.?\d*)Mb", content)
    metrics['quality_cutoffs'] = {}
    for cutoff in quality_cutoffs:
        quality_score, count, percentage, mb = cutoff
        metrics['quality_cutoffs'][f'Q{quality_score}'] = {
            'count': int(count),
            'percentage': float(percentage),
            'megabases': float(mb.replace('Mb', ''))
        }

    # Parsing top reads information (mean basecall quality)
    top_reads = re.findall(r"(\d+):\s*(\d+\.\d+)\s*\((\d+)\)", content)
    metrics['top_reads_quality'] = [(int(top[0]), float(top[1]), int(top[2])) for top in top_reads]

    return metrics


def generate_qc_report(metrics, api_key, species="Human", tissue="Skin", model="deepseek-reasoner", output_text_file='qc_report.txt'):
    """
    Generate a quality control report from the parsed NanoStats data.
    
    Parameters:
    metrics (dict): Dictionary containing parsed sequencing metrics.
    api_key (str): API key for the model.
    species (str): The species of the sample (default is "Human").
    tissue (str): The tissue type of the sample (default is "Skin").
    model (str): The model to use (default is "deepseek-reasoner").
    output_text_file (str): Path to save the generated QC report.
    
    Returns:
    str: The content of the generated QC report.
    """
    client = OpenAI(api_key=api_key, base_url='https://api.deepseek.com')

    # Construct the system prompt based on the parsed data and provided species/tissue
    system_prompt = (
        f"Generate a quality control report for the following single-cell Nanopore long-read RNA-seq data from {species} {tissue} tissue. "
        "The report should evaluate each of the following metrics and provide a conclusion about the data quality, "
        "with suggestions for improving the quality if necessary:\n\n"
        f"Mean Read Length: {metrics['mean_read_length']} bp\n"
        f"Mean Read Quality: {metrics['mean_read_quality']}\n"
        f"Median Read Length: {metrics['median_read_length']} bp\n"
        f"Median Read Quality: {metrics['median_read_quality']}\n"
        f"Number of Reads: {metrics['number_of_reads']}\n"
        f"Read Length N50: {metrics['read_length_n50']} bp\n"
        f"STDEV Read Length: {metrics['stdev_read_length']} bp\n"
        f"Total Bases: {metrics['total_bases']} bases\n\n"
        "Quality Cutoffs:\n"
    )

    for cutoff, data in metrics['quality_cutoffs'].items():
        system_prompt += f"> {cutoff}: {data['count']} reads ({data['percentage']}%) {data['megabases']}Mb\n"

    system_prompt += "\nTop 5 Highest Mean Basecall Quality Scores and Their Read Lengths:\n"
    for i, (read_num, quality, length) in enumerate(metrics['top_reads_quality']):
        system_prompt += f"{i + 1}: {read_num} (Quality: {quality}, Length: {length} bp)\n"

    # Sending the request to DeepSeek API
    messages = [{"role": "system", "content": system_prompt}]
    
    completion = client.chat.completions.create(
        model=model,
        messages=messages
    )

    # Get the response content from the model
    response_content = completion.choices[0].message.content
    
    # Save the result to the output text file
    with open(output_text_file, 'w') as f:
        f.write(response_content)
    
    print(f"Quality control report saved to {output_text_file}")

    return response_content


def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Generate a QC report for sequencing data.")
    
    # Arguments
    parser.add_argument('--file_path', required=True, help='Path to the NanoStats text file')
    parser.add_argument('--api_key', required=True, help='API key for DeepSeek model')
    parser.add_argument('--output_text_file', required=True, help='Path to save the output QC report')
    parser.add_argument('--species', default="Human", help='The species of the sample (default: Human)')
    parser.add_argument('--tissue', default="Skin", help='The tissue type of the sample (default: Skin)')
    parser.add_argument('--model', default="deepseek-reasoner", help='The model to use (default: deepseek-reasoner)')

    # Parse the arguments
    args = parser.parse_args()

    # Parse the NanoStats file and extract metrics
    metrics = parse_nanostats_file(args.file_path)

    # Generate the QC report
    generate_qc_report(
        metrics, 
        api_key=args.api_key, 
        species=args.species, 
        tissue=args.tissue, 
        model=args.model, 
        output_text_file=args.output_text_file
    )


if __name__ == '__main__':
    main()

#python /data/haowu/sc_long_v2/pipline_v2/bin/read_qc_deekseep.py --file_path /data/haowu/sc_long_v2/pipline_v2/demo/all_samples_NanoStats.txt \
                          # --api_key sk-403bc4c7826a475e9646978562cb324a \
                          # --output_text_file /data/haowu/sc_long_v2/pipline_v2/demo/qc_report.txt \
                          # --species "Human" --tissue "Skin"
