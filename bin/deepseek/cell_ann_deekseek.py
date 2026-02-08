import pandas as pd
from openai import OpenAI
import argparse

prompt_basic = """ Annotate cell clusters for {species_type} {tissue_type} restricting output STRICTLY to Major Cell Lineages. PROHIBITED: Subtypes, subsets, activation states, anatomical location specificity, or disease states. REQUIRED: Output must map to the highest-level biological category applicable to {tissue_type}.
Each cluster and its corresponding cell type should be listed on a separate line in the format: Cluster X: [Cell Type].
Single line only give one result!
If markers are unrecognized, you MUST predict the closest Major Cell Lineage based on functional enrichment or gene similarity. The output 'Unknown' or 'Insufficient_Data' is STRICTLY PROHIBITED.
Perform the following cell type annotation process on the provided gene list for each cluster. Conduct this process 5 times in parallel,
 independently. After the 5 runs, determine the final annotation for each cluster by selecting the most frequent annotation. If there's a tie in the frequency of annotations,
  choose the annotation with the highest average confidence level among the tied annotations. If both frequency and confidence levels are tied, 
  prioritize the annotation based on this confidence level order: [3] > [2] > [1] > [0]. For each cluster (represented by 'X'), 
  output only the final annotation in the format 'ClusterX:Cell_TYPE[i]', where Cell_TYPE is the determined cell type and 'i' is the confidence level. 
  Do not include any additional text, formatting, or explanations. The output must strictly adhere to this format to be correctly processed.
Provide your response without bullet points, indentation, or other formatting.
1.  Initialization
    1.1  Calculate Genes_Per_Cluster: Count genes provided per cluster.
    1.2  Determine Top_N Selection:
        (1) Genes_Per_Cluster < 5: Top_N = All Genes; If no highly specific core markers, proceed immediately to Functional Enrichment and Gene Similarity analysis to force a prediction (see Section 2.4).
        (2) 5 ≤ Genes_Per_Cluster < 10: Top_N = 5; Apply stricter confidence criteria (see Section 5).
        (3) 10 ≤ Genes_Per_Cluster < 15: Top_N = 10.
        (4) Genes_Per_Cluster ≥ 15: Top_N = 15; Consider Top_30 for exclusion, multi-lineage, disease, tissue-specific, closely related, and progenitor markers.
2.  Output Specification
    2.1  Format: ClusterX: Primary_Cell_Type [Confidence_Level]
    2.2  Examples:
        Cluster2: Hepatocyte[3]
        Cluster5: Endothelial/Fibroblast[2]
        Cluster8: Unknown_Immune*[1]
    2.3  Multi-Lineage Format: CellTypeA/CellTypeB[Confidence] if evidence score supports multiple possible annotations.
    2.4 Forced Prediction Protocol (No Unknowns):
If core markers are absent or unrecognizable, apply the following steps sequentially to assign a Major Cell Lineage:
(1) Functional Enrichment: Predict lineage based on gene functions/pathways (e.g., GO terms implied by the list).
(2) Gene Similarity: If (1) fails, evaluate gene list similarity to known lineage profiles.
(3) Output Requirement: Output the closest matching Major Cell Lineage with confidence level [0].
(4) Prohibition: NEVER output 'Unknown', 'Unassigned', or 'Insufficient_Data'.
3.  Primary cell type assignment (Cell Ontology - CL)
    3.1  Specificity Constraints:
        (1) Maximum: CL level 3 (Exceptions: disease, tissue-specific, progenitors. See Sections 7-9).
        (2) Minimum: CL level 2. If Genes_Per_Cluster < 10, relax to CL level 1 only if no CL level 2 terms are supported.
    3.2  Core Marker Validation:
        (1) Positive (Required):
            a) Genes_Per_Cluster < 10: ≥ 2 highly specific core markers in Top_N. If only 1, reduce confidence by 1 level.
            b) Genes_Per_Cluster ≥ 10: ≥ 3 core markers in Top_N. If exactly 3, append *.
        (2) Negative (Exclusion Markers):
            a) Genes_Per_Cluster < 15: No exclusion markers in the entire gene list.
            b) Genes_Per_Cluster ≥ 15: Max 1 exclusion marker in Top_30 only if: (A) Not in authoritative exclusion lists, AND (B) Rank > Top_10.
            c) Exclusion Marker Penalty: If found in Top_N, reduce confidence or append *.
4.  Multi-Lineage Resolution
    4.1  Trigger multi-lineage annotation when:
        (1) Top_N contains ≥1 definitive marker from ≥2 distinct lineages.
        (2) No single lineage surpasses 60% evidence score.
    4.2  Evidence Score Calculation:
        (1) Evidence Score = Σ (marker_rank_weight × specificity_score)
        (2) Rank Weights (based on gene rank in Top_N or Top_30): Ranks 1-5: 3×, Ranks 6-10: 2×, Ranks 11-30: 1×.
        (3) Specificity Scores: Core markers: 3, Tissue-specific marker: 2, Pan-lineage marker: 1, Disease-specific, progenitor, or closely related: 2.
        (4) If the highest scoring lineage has < 60% of total possible score, use slash notation to reflect ambiguity.
5.  Confidence Assignment
    5.1  [3] (Unambiguous): Unique lineage markers present.
    5.2  [2] (Probable): Supporting markers present, no conflicts.
    5.3  [1] (Tentative): Weak markers, no exclusions.
    5.4 [0] (Forced Prediction): Annotation derived solely from Functional Enrichment or Gene Similarity due to lack of specific core markers.
    5.5  [*] (Requires Validation): Append * to any annotation triggering a validation flag.
6.  Validation Checks
    6.1  Orthology Verification (Cross-species):
        (1) Forward and reverse orthologs must map to the same cell type.
        (2) Use Ensembl Compara or OMA Browser with ≥ 90% orthology confidence.
        (3) Invalid or ambiguous orthology: Map to the broad lineage class (e.g., Animal_Cell) or the closest predicted lineage based on gene function with confidence [0]*.
    6.2  Tissue Plausibility Check:
        (1) Annotated cell type must be biologically plausible in {tissue_type}.
        (2) For rare or novel tissues, allow relaxed plausibility but append *.
        (3) Inconsistent annotations (e.g,, neurons in liver) must be flagged with * or rejected.
    6.3  Mitochondrial Gene Inclusion:
        (1) MT genes allowed only if:
            a) At least 2 structural markers consistent with high-mito content cells are present (e.g., FOXJ1+ for ciliated cells).
            b) MT% < threshold (Animals: 15%).
7.  Special cases
    7.1  Disease-specific annotation (including cancer):
        (1) Indicators: Disease-associated markers (e.g., oncogenes, fibrosis markers, inflammatory cytokines), aberrant gene expression, or lineage reprogramming signals.
        (2) Validation Requirements:
            a) At least 1 marker must be present in a validated disease gene database (e.g., COSMIC for cancer, DisGeNET for other diseases).
            b) At least 2 supporting markers should be associated with the same disease or condition in literature or knowledge bases (e.g., CancerSEA, Harmonizome, HuDiNe).
            c) Disease-specific annotations must also pass tissue plausibility and exclusion checks.
        (3) Format: Tissue_DiseaseSubtype[Confidence]* (e.g., Lung_Adenocarcinoma[2]*, Liver_FibrosisAssociated[1]*, Intestinal_InflammatoryCell[2]*)
        (4) Note: Use the most specific disease label possible, but prefer general disease categories if gene evidence is weak or ambiguous.
    7.2  Tissue-specific annotation:
        (1) Markers must be highly enriched in {tissue_type} and supported by ≥ 2 references (e.g., PanglaoDB, Human Protein Atlas).
        (2) Format: TissueSpecific_CellType[Confidence] or {tissue_type}_Specific_Cell[Confidence]*
    7.3  Progenitor cell annotation:
        (1) Indicators: stemness signatures, mild cell cycle markers, lineage commitment genes.
        (2) Avoid over-specificity (e.g., use "Myeloid_Progenitor" over "Common_Monocyte_Progenitor").
        (3) Format: Lineage_Progenitor[Confidence]* or Progenitor_Cell[1]*
8.  Automatic Rejection Conditions
    Reject annotation automatically if:
    (1) Any sub-type marker in Top_20.
    (2) Cell cycle genes dominate Top_N.
    (3) 30% ribosomal genes without clear secretory function.
(4) Core markers conflict across annotations.
(5) The output contains the string "Unknown", "Unassigned", or "Insufficient". (Must re-run prediction via Section 2.4).
9.  Knowledge-base Prioritization
    9.1  Highest: Perfect CL match.
    9.2  Second: ≥ 3 publications support marker combination.
    9.3  Third: Pan-tissue/low-specificity annotations if necessary.
    9.4  Blacklist: Maintain a list of commonly misannotated terms (e.g., "Fibroblast" in brain). Blacklisted terms require ≥ 1 tissue-specific marker in Top_5 to be accepted. ** Maintain a list of commonly misannotated terms requiring additional validation.
"""

def cell_ann_gpt(args):
    """
    Automatically annotate cell types using DeepSeek API based on gene expression data.
    
    Parameters:
    args (Namespace): Arguments passed from the command line.
    
    Returns:
    str: Cell type annotation results returned by the DeepSeek API.
    list: Conversation history containing requests and responses.
    """
    client = OpenAI(
        api_key=args.api_key,
        base_url='https://api.deepseek.com'
    )

    # Read CSV file
    df = pd.read_csv(args.csv_file_path, header=None)

    pre_prompt = (f"Annotate cell clusters for {args.species_type} {args.tissue_type} restricting output STRICTLY to Major Cell Lineages. "
                  f"PROHIBITED: Subtypes, subsets, activation states, anatomical location specificity, or disease states. "
                  f"REQUIRED: Output must map to the highest-level biological category applicable to {args.tissue_type}."
                  f"Each cluster and its corresponding cell type should be listed on a separate line in the format: Cluster X: [Cell Type]."
                  f"Single line only give one result! \n")
    
    # Conversation history
    dialogue_history = []

    # Extract cluster and gene information
    clusters = df.iloc[:, 0].tolist()
    genes = df.iloc[:, 1].str.split(',').tolist()
    
    # Build system prompt
    system_prompt = pre_prompt + prompt_basic
    messages = [{"role": "system", "content": system_prompt}]  # Initialize messages list

    user_prompt_content = ""  # Accumulate prompts for all clusters
    for cluster, gene_list in zip(clusters, genes):
        # Limit the number of genes to display
        gene_list = gene_list[:args.max_genes_to_display]

        # Build prompt for current cluster
        cluster_prompt = f"Cluster {cluster}: {', '.join(gene_list)}\n"
        user_prompt_content += cluster_prompt

    messages.append({"role": "user", "content": user_prompt_content})  # Add user content

    # Send the final prompt to DeepSeek API
    completion = client.chat.completions.create(
        model=args.model,
        messages=messages
    )

    # Record API response
    response_content = completion.choices[0].message.content
    dialogue_history.append({"prompt": messages, "response": response_content, "role": "assistant"})  # Record the conversation

    # Parse the response content into annotations for each cluster
    annotated_cells = []
    lines = response_content.strip().split('\n')
    cluster_count = 0

    for line in lines:
        if line.startswith("Cluster"):
            try:
                parts = line.split(":")
                cluster_id = parts[0].strip()
                cell_type = parts[1].strip()
                annotated_cells.append({"Cluster": cluster_id, "Cell Type": cell_type})
            except IndexError:
                # If cannot split into two columns, automatically generate a Cluster column
                cluster_id = f"Cluster{cluster_count}"
                cell_type = line.strip()
                cluster_count += 1
                annotated_cells.append({"Cluster": cluster_id, "Cell Type": cell_type})

    # Save results to CSV file
    annotated_df = pd.DataFrame(annotated_cells)
    annotated_df.to_csv(args.output_csv, index=False)
    print(f"Annotation results saved to {args.output_csv}")

    # Return DeepSeek API response and conversation history
    return response_content, dialogue_history


def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Automatically annotate cell types using DeepSeek API")
    
    # Required arguments
    parser.add_argument('--csv_file_path', required=True, help='Path to CSV file containing cluster and gene information')
    parser.add_argument('--api_key', required=True, help='DeepSeek API key')
    parser.add_argument('--output_csv', required=True, help='Output CSV file path')
    
    # Optional arguments (with defaults)
    parser.add_argument('--tissue_type', default='Skin', help='Tissue type, default is "Skin"')
    parser.add_argument('--species_type', default='human', help='Species type, default is "human"')
    parser.add_argument('--max_genes_to_display', type=int, default=30, 
                       help='Maximum number of genes to display per cluster, default is 30')
    parser.add_argument('--model', default='deepseek-reasoner', 
                       help='DeepSeek model to use, default is "deepseek-reasoner"')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Call cell annotation function
    annotation_result, dialogue_history = cell_ann_gpt(args)
    
    return annotation_result


if __name__ == '__main__':
    main()
