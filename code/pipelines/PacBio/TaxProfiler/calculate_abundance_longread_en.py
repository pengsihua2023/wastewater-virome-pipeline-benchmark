#!/usr/bin/env python3
"""
Long-read relative abundance calculation script (based on the Kraken2 report)

For long-read data (Nanopore/PacBio), Kraken2 classification is already accurate enough,
No Bracken statistical correction is needed. Relative abundance (RPM) is extracted directly from the Kraken2 report.

Usage:
    python calculate_abundance_longread_en.py \
        --kraken results/kraken2/sample.report \
        --output sample_abundance.tsv
"""

import argparse
import pandas as pd
import sys
from pathlib import Path


def parse_kraken_report(kraken_file, taxonomy_level='S'):
    """
    Parse Kraken2 report file, extract abundance at specified taxonomic level
    
    Kraken2 report format:
    percentage  clade_reads  taxon_reads  rank  taxid  name
    """
    print(f"📖 Reading Kraken2 report: {kraken_file}")
    
    try:
        # Read Kraken2 report (no header, tab-separated)
        df = pd.read_csv(
            kraken_file, 
            sep='\t', 
            header=None,
            names=['percentage', 'clade_reads', 'taxon_reads', 'rank', 'taxid', 'name']
        )
        
        # Clean species names (remove leading spaces and rank markers)
        df['name'] = df['name'].str.strip()
        
        # Filter to specified taxonomic level
        # S = Species, G = Genus, F = Family, O = Order, C = Class, P = Phylum, D = Domain
        df_filtered = df[df['rank'] == taxonomy_level].copy()
        
        print(f"✅ Found {len(df_filtered)} {get_level_name(taxonomy_level)}-level taxonomic units")
        
        return df_filtered
        
    except Exception as e:
        print(f"❌ Error: Unable to read Kraken2 report - {e}")
        sys.exit(1)


def get_level_name(level_code):
    """Get taxonomic level name"""
    level_map = {
        'D': 'domain',
        'P': 'phylum',
        'C': 'class',
        'O': 'order',
        'F': 'family',
        'G': 'genus',
        'S': 'species'
    }
    return level_map.get(level_code, level_code)


def calculate_abundance(kraken_df):
    """
    Compute the relative abundance metric RPM (Reads Per Million)
    
    For long-read data:
    - Use clade_reads from the Kraken2 report (reads including child taxa)
    - percentage is already a relative abundance (0-100)
    """
    print("\n📊 Calculating relative abundance (RPM)...")
    
    # Compute the total number of reads (using the root node reads)
    # Note: the percentages in the Kraken2 report are relative to the total reads
    # We need to back-calculate the total number of reads from the percentage
    if len(kraken_df) > 0:
        # Use clade_reads and percentage to compute the total reads
        max_clade = kraken_df['clade_reads'].max()
        max_pct = kraken_df[kraken_df['clade_reads'] == max_clade]['percentage'].values[0]
        if max_pct > 0:
            total_reads = int(max_clade * 100 / max_pct)
        else:
            total_reads = kraken_df['clade_reads'].sum()
    else:
        print("⚠️  Warning: No classification data found")
        return pd.DataFrame()
    
    print(f"   Total reads (estimated): {total_reads:,}")
    
    # Create the result list
    results = []
    
    for idx, row in kraken_df.iterrows():
        species_name = row['name']
        taxonomy_id = row['taxid']
        assigned_reads = row['clade_reads']  # Use clade_reads
        fraction = row['percentage'] / 100  # convert to 0-1
        
        # Compute RPM (relative abundance)
        rpm = (assigned_reads / total_reads) * 1_000_000 if total_reads > 0 else 0
        
        results.append({
            'Species': species_name,
            'Taxonomy_ID': taxonomy_id,
            'Assigned_Reads': int(assigned_reads),
            'Fraction': fraction,
            'RPM': rpm
        })
    
    result_df = pd.DataFrame(results)
    
    # Sort by RPM in descending order
    result_df = result_df.sort_values('RPM', ascending=False)
    
    print(f"✅ Calculation complete, total {len(result_df)} species")
    
    return result_df


def main():
    parser = argparse.ArgumentParser(
        description='Long-read data abundance calculation (based on Kraken2 report)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # Single sample
  python calculate_abundance_longread_en.py \\
    --kraken results/kraken2/sample1.report \\
    --output sample1_abundance.tsv
  
  # Specify taxonomic level
  python calculate_abundance_longread_en.py \\
    --kraken results/kraken2/sample1.report \\
    --level G \\
    --output sample1_abundance_genus.tsv

Taxonomic level codes:
  D = Domain
  P = Phylum
  C = Class
  O = Order
  F = Family
  G = Genus
  S = Species (default)
        """
    )
    
    parser.add_argument('--kraken', required=True,
                       help='Kraken2 report file path (*.report)')
    parser.add_argument('--output', required=True,
                       help='Output file path')
    parser.add_argument('--level', default='S',
                       help='Taxonomic level (default: S=species)')
    args = parser.parse_args()
    
    print("="*60)
    print("🧬 Long-read Viral Abundance Calculator (RPM)")
    print("   based on the Kraken2 report (no Bracken)")
    print("="*60)
    
    # Read the Kraken2 report
    kraken_df = parse_kraken_report(args.kraken, args.level)
    
    if len(kraken_df) == 0:
        print("\n⚠️  Warning: No classification data found")
        # Create an empty result file
        pd.DataFrame(columns=['Species', 'Taxonomy_ID', 'Assigned_Reads', 
                             'Fraction', 'RPM']).to_csv(
            args.output, sep='\t', index=False
        )
        sys.exit(0)
    
    # Compute the relative abundance
    result_df = calculate_abundance(kraken_df)
    
    # Save results
    print(f"\n💾 Saving results to: {args.output}")
    result_df.to_csv(args.output, sep='\t', index=False)
    
    # Print statistical summary
    print("\n" + "="*60)
    print("📈 Statistical Summary")
    print("="*60)
    print(f"Detected {get_level_name(args.level)} count: {len(result_df)}")
    print(f"Total reads: {result_df['Assigned_Reads'].sum():,}")
    print(f"\nTop 10 most abundant {get_level_name(args.level)} (by RPM):")
    print("-"*60)
    
    top10 = result_df.head(10)
    for idx, row in top10.iterrows():
        print(f"{row['Species'][:50]:50s} RPM: {row['RPM']:10.2f}")
    
    print("\n✅ Complete!")


if __name__ == '__main__':
    main()

