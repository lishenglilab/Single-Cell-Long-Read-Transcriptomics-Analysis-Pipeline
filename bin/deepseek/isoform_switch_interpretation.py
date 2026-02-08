#!/usr/bin/env python3
"""
Isoform Switch Biological Interpretation using DeepSeek API
Analyzes isoform switching consequences and provides biological interpretation
"""

import pandas as pd
import argparse
import sys
import os
from openai import OpenAI
import json
from datetime import datetime

def parse_consequences_file(consequences_file):
    """
    Parse switch_consequences.txt file and extract significant consequences
    
    Parameters:
    consequences_file (str): Path to switch_consequences.txt file
    
    Returns:
    dict: Parsed consequences data and gene information
    """
    try:
        df = pd.read_csv(consequences_file, sep="\t")
        
        # Filter rows with isoformsDifferent = True
        df_true = df[df['isoformsDifferent'] == True]
        
        if df_true.empty:
            return {
                'status': 'no_significant',
                'message': 'No rows with isoformsDifferent = True found.'
            }
        
        # Get first gene information
        first_gene_row = df_true.iloc[0]
        gene_name = first_gene_row['gene_name']
        condition_1 = first_gene_row['condition_1']
        condition_2 = first_gene_row['condition_2']
        isoform_up = first_gene_row['isoformUpregulated']
        isoform_down = first_gene_row['isoformDownregulated']
        
        # Get all true consequences for this gene
        gene_consequences = df_true[df_true['gene_name'] == gene_name]
        true_consequences = []
        
        for _, row in gene_consequences.iterrows():
            if row['isoformsDifferent'] == True and pd.notna(row['switchConsequence']):
                feature = row['featureCompared']
                consequence = row['switchConsequence']
                true_consequences.append({
                    'feature': feature,
                    'consequence': consequence
                })
        
        return {
            'status': 'success',
            'gene_info': {
                'gene_name': gene_name,
                'condition_1': condition_1,
                'condition_2': condition_2,
                'isoform_up': isoform_up,
                'isoform_down': isoform_down
            },
            'consequences': true_consequences,
            'total_consequences': len(true_consequences)
        }
        
    except Exception as e:
        return {
            'status': 'error',
            'message': f'Error parsing consequences file: {str(e)}'
        }

def generate_prompt(gene_data, study_description):
    """
    Generate prompt for DeepSeek API
    
    Parameters:
    gene_data (dict): Gene information and consequences
    study_description (str): Study context description
    
    Returns:
    str: Formatted prompt for API
    """
    gene_info = gene_data['gene_info']
    consequences = gene_data['consequences']
    
    # Format consequences list
    consequences_list = "\n".join([
        f"  • {c['feature']}: {c['consequence']}" 
        for c in consequences
    ])
    
    prompt = f"""You are a senior bioinformatics expert and molecular biologist. Please provide a professional, concise biological interpretation based on the following research context and experimental results.

Research Context:
{study_description}

Experimental Results:
Gene {gene_info['gene_name']} shows isoform switching between {gene_info['condition_1']} and {gene_info['condition_2']} conditions, with transcript {gene_info['isoform_up']} upregulated in {gene_info['condition_2']} and transcript {gene_info['isoform_down']} dominant in {gene_info['condition_1']}.

The detected functional consequences of this isoform switch include:
{consequences_list}

Please provide interpretation covering the following aspects:
1. Known function of this gene (if known)
2. How these consequences (e.g., ORF shortening, domain loss, IDR loss) might affect protein function
3. Potential implications for cellular function in the disease context (AD, atopic dermatitis)
4. Possible regulatory mechanisms (e.g., splicing factor changes)

Please respond in English and ensure the explanation is concise, professional, and suitable for biologists.
"""
    return prompt

def call_deepseek_api(prompt, api_key, model="deepseek-chat"):
    """
    Call DeepSeek API for biological interpretation
    
    Parameters:
    prompt (str): Prompt for the model
    api_key (str): DeepSeek API key
    model (str): Model to use
    
    Returns:
    dict: API response and status
    """
    try:
        client = OpenAI(api_key=api_key, base_url='https://api.deepseek.com')
        
        response = client.chat.completions.create(
            model=model,
            messages=[
                {
                    "role": "system", 
                    "content": "You are a senior bioinformatics expert and molecular biologist specializing in transcriptomics and disease mechanisms."
                },
                {
                    "role": "user", 
                    "content": prompt
                }
            ],
            temperature=0.7,
            max_tokens=2000
        )
        
        interpretation = response.choices[0].message.content
        
        return {
            'status': 'success',
            'interpretation': interpretation,
            'model': model,
            'tokens_used': response.usage.total_tokens if hasattr(response, 'usage') else 'unknown'
        }
        
    except Exception as e:
        return {
            'status': 'error',
            'message': f'API call failed: {str(e)}'
        }

def main():
    parser = argparse.ArgumentParser(
        description='Generate biological interpretation for isoform switching using DeepSeek API'
    )
    
    # Required arguments
    parser.add_argument('--consequences_file', required=True, 
                       help='Path to switch_consequences.txt file')
    parser.add_argument('--output_report', required=True,
                       help='Path to save the biological interpretation report')
    parser.add_argument('--output_api', required=True,
                       help='Path to save the API response details')
    
    # Optional arguments
    parser.add_argument('--api_key', default='',
                       help='DeepSeek API key (required for API-based interpretation)')
    parser.add_argument('--study_description', default='',
                       help='Study context description')
    parser.add_argument('--model', default='deepseek-chat',
                       help='DeepSeek model to use (default: deepseek-chat)')
    parser.add_argument('--save_prompt', action='store_true',
                       help='Save the generated prompt to a file')
    
    args = parser.parse_args()
    
    # Check if API key is provided
    if not args.api_key or args.api_key == '':
        print("ERROR: DeepSeek API key is not provided. Please set --api_key parameter.")
        with open(args.output_report, 'w') as f:
            f.write("Biological interpretation skipped: API key not provided.\n")
            f.write("To generate biological interpretation, please provide a DeepSeek API key.\n")
        with open(args.output_api, 'w') as f:
            f.write(json.dumps({
                'status': 'skipped',
                'message': 'API key not provided',
                'timestamp': datetime.now().isoformat()
            }, indent=2))
        sys.exit(0)
    
    # Set default study description if not provided
    if not args.study_description:
        args.study_description = """This is a single-cell full-length transcriptome sequencing analysis comparing Disease and Normal conditions. 
The study aims to identify differentially expressed transcripts and alternative splicing events in Disease, 
exploring disease-associated isoform switching and its functional consequences."""
    
    print(f"Parsing consequences file: {args.consequences_file}")
    
    # Parse the consequences file
    gene_data = parse_consequences_file(args.consequences_file)
    
    if gene_data['status'] == 'error':
        print(f"ERROR: {gene_data['message']}")
        with open(args.output_report, 'w') as f:
            f.write(f"Error parsing consequences file: {gene_data['message']}")
        sys.exit(1)
    
    if gene_data['status'] == 'no_significant':
        print("No significant isoform switch consequences found.")
        with open(args.output_report, 'w') as f:
            f.write("No significant isoform switch consequences (isoformsDifferent = True) were found to interpret.")
        with open(args.output_api, 'w') as f:
            f.write(json.dumps({
                'status': 'no_significant',
                'message': 'No significant consequences found',
                'timestamp': datetime.now().isoformat()
            }, indent=2))
        sys.exit(0)
    
    # Generate prompt
    prompt = generate_prompt(gene_data, args.study_description)
    
    # Save prompt if requested
    if args.save_prompt:
        prompt_file = args.output_report.replace('.txt', '_prompt.txt')
        with open(prompt_file, 'w') as f:
            f.write(prompt)
        print(f"Prompt saved to: {prompt_file}")
    
    # Call DeepSeek API
    print("Calling DeepSeek API...")
    api_result = call_deepseek_api(prompt, args.api_key, args.model)
    
    if api_result['status'] == 'success':
        interpretation = api_result['interpretation']
        
        # Save interpretation report
        with open(args.output_report, 'w') as f:
            f.write(interpretation)
        
        # Save API response details
        api_details = {
            'timestamp': datetime.now().isoformat(),
            'status': 'success',
            'gene_info': gene_data['gene_info'],
            'total_consequences': gene_data['total_consequences'],
            'model': api_result['model'],
            'tokens_used': api_result['tokens_used'],
            'prompt': prompt,
            'response': interpretation
        }
        
        with open(args.output_api, 'w') as f:
            json.dump(api_details, f, indent=2)
        
        print(f"Biological interpretation generated for gene {gene_data['gene_info']['gene_name']}")
        print(f"Report saved to: {args.output_report}")
        
    else:
        print(f"API call failed: {api_result['message']}")
        
        with open(args.output_report, 'w') as f:
            f.write(f"API call failed: {api_result['message']}\n")
            f.write("\nNo biological interpretation was generated due to API error.\n")
        
        with open(args.output_api, 'w') as f:
            f.write(json.dumps({
                'status': 'api_error',
                'message': api_result['message'],
                'timestamp': datetime.now().isoformat()
            }, indent=2))
    
    # Create gene summary file
    summary_file = 'gene_summary.txt' 
    summary = f"""Gene Information Summary:
    - Gene Name: {gene_data['gene_info']['gene_name']}
    - Comparison: {gene_data['gene_info']['condition_1']} vs {gene_data['gene_info']['condition_2']}
    - Upregulated Transcript: {gene_data['gene_info']['isoform_up']}
    - Downregulated Transcript: {gene_data['gene_info']['isoform_down']}
    - Significant Consequences Detected: {gene_data['total_consequences']}
    - Analysis Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
    """
    
    with open(summary_file, 'w') as f:
        f.write(summary)
    
    print(f"Analysis completed. Summary saved to: {summary_file}")

if __name__ == '__main__':
    main()
