#!/usr/bin/env python3
"""
Compute the viral relative abundance metric: RPM (Reads Per Million)

RPM (Reads Per Million): number of reads assigned to the species per million reads
This is a relative abundance metric used to compare the relative abundance of the same virus across samples

Usage:
    python calculate_abundance_en.py --bracken results/bracken/sample_bracken.tsv \
                                      --kraken results/kraken2/sample.report \
                                      --output results/abundance_RPM_RPKM.tsv
"""

import argparse
import pandas as pd
import sys
from pathlib import Path


def parse_bracken_output(bracken_file):
    """
    Parse Bracken output file
    
    Bracken output format:
    name    taxonomy_id    taxonomy_lvl    kraken_assigned_reads    added_reads    new_est_reads    fraction_total_reads
    """
    print(f"📖 Reading Bracken file: {bracken_file}")
    
    try:
        df = pd.read_csv(bracken_file, sep='\t')
        print(f"✅ Successfully read {len(df)} taxonomic units")
        
        # Show the column names for debugging
        print(f"   Columns found: {list(df.columns)}")
        
        # Handle different column name variants
        # Bracken files may have different column name formats
        column_mapping = {
            'name': ['name', 'Name', 'species', 'Species'],
            'taxonomy_id': ['taxonomy_id', 'taxonomyId', 'taxid', 'Taxonomy_ID'],
            'new_est_reads': ['new_est_reads', 'newEstReads', 'reads', 'Reads'],
            'fraction_total_reads': ['fraction_total_reads', 'fractionTotalReads', 'fraction', 'Fraction']
        }
        
        # Standardize the column names
        for standard_name, possible_names in column_mapping.items():
            for possible_name in possible_names:
                if possible_name in df.columns:
                    if standard_name != possible_name:
                        df[standard_name] = df[possible_name]
                    break
        
        # Verify that the required columns exist
        required_cols = ['name', 'taxonomy_id', 'new_est_reads', 'fraction_total_reads']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            print(f"⚠️  Warning: Missing columns: {missing_cols}")
            print(f"   Available columns: {list(df.columns)}")
            print(f"   First few rows:")
            print(df.head())
        
        return df
    except Exception as e:
        print(f"❌ Error: Unable to read Bracken file - {e}")
        sys.exit(1)


def parse_kraken_report(kraken_file):
    """
    Parse Kraken2 report file for genome length information
    
    Kraken2 report format (tab-separated):
    percentage    num_covered    num_assigned    rank    taxid    name
    """
    print(f"📖 Reading Kraken2 report: {kraken_file}")
    
    try:
        # Kraken2 report has no header
        # Try different numbers of columns (some reports may have extra columns)
        df = pd.read_csv(kraken_file, sep='\t', header=None)
        
        # Determine the format from the number of columns
        if len(df.columns) >= 6:
            df.columns = ['percentage', 'num_covered', 'num_assigned', 
                         'rank', 'taxid', 'name'] + [f'extra_{i}' for i in range(6, len(df.columns))]
        elif len(df.columns) == 5:
            df.columns = ['percentage', 'num_covered', 'num_assigned', 'rank', 'taxid']
            df['name'] = ''  # if there is no name column, set it to empty
        else:
            print(f"⚠️  Warning: Unexpected number of columns ({len(df.columns)}) in Kraken2 report")
            print(f"   First row: {df.iloc[0].tolist()}")
        
        # Clean species names (remove leading spaces and rank markers)
        if 'name' in df.columns:
            df['name'] = df['name'].str.strip()
        
        print(f"✅ Successfully read Kraken2 report with {len(df)} entries")
        return df
    except Exception as e:
        print(f"❌ Error: Unable to read Kraken2 report - {e}")
        sys.exit(1)


def calculate_rpm(bracken_df, kraken_df):
    """
    Compute the relative abundance metric RPM (Reads Per Million)
    
    RPM = (assigned_reads / total_reads) × 1,000,000
    
    This is a relative abundance metric used to compare the relative abundance of the same virus across samples
    """
    print("\n📊 Calculating relative abundance (RPM)...")
    
    # Verify that the required columns exist
    if 'new_est_reads' not in bracken_df.columns:
        print(f"❌ Error: Column 'new_est_reads' not found in Bracken data")
        print(f"   Available columns: {list(bracken_df.columns)}")
        sys.exit(1)
    
    # Compute the total number of reads
    total_reads = bracken_df['new_est_reads'].sum()
    
    print(f"   Total reads: {total_reads:,}")
    
    # Verify the data
    if total_reads == 0:
        print("⚠️  Warning: the total number of reads is 0, please check the input file")
        sys.exit(1)
    
    # Verify the sum of fraction_total_reads (should be close to 1.0)
    if 'fraction_total_reads' in bracken_df.columns:
        fraction_sum = bracken_df['fraction_total_reads'].sum()
        print(f"   Fraction sum: {fraction_sum:.6f} (should be close to 1.0)")
        if abs(fraction_sum - 1.0) > 0.1:
            print(f"   ⚠️  Warning: the Fraction sum deviates substantially from 1.0, there may be a problem")
    
    # Create the result list
    results = []
    
    # Statistics
    processed_count = 0
    zero_reads_count = 0
    
    for idx, row in bracken_df.iterrows():
        # Safely get the column value; use the default if it does not exist
        species_name = row.get('name', f'Unknown_{idx}')
        taxonomy_id = row.get('taxonomy_id', 'NA')
        assigned_reads = row['new_est_reads']  # this column must exist
        
        # Skip species with 0 reads (optional)
        if assigned_reads == 0:
            zero_reads_count += 1
            continue
        
        # Use the fraction provided by Bracken; compute it if absent
        if 'fraction_total_reads' in bracken_df.columns:
            fraction = row.get('fraction_total_reads', assigned_reads / total_reads if total_reads > 0 else 0)
        else:
            fraction = assigned_reads / total_reads if total_reads > 0 else 0
        
        # Debug information (first few only)
        if processed_count < 3:
            print(f"   Processing {processed_count+1}: {species_name[:40]:40s}, reads: {assigned_reads:>12,}, fraction: {fraction:.6f}")
        
        # Compute RPM (relative abundance)
        rpm = (assigned_reads / total_reads) * 1_000_000 if total_reads > 0 else 0
        
        processed_count += 1
        
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
    
    # Verify the results
    result_total_reads = result_df['Assigned_Reads'].sum()
    result_fraction_sum = result_df['Fraction'].sum()
    
    print(f"✅ Calculation complete, total {len(result_df)} species")
    print(f"   Result check:")
    print(f"   - Species processed: {processed_count}")
    print(f"   - Skipped species with 0 reads: {zero_reads_count}")
    print(f"   - Total reads in the result reads: {result_total_reads:,}")
    print(f"   - Sum of Fraction in the result: {result_fraction_sum:.6f}")
    
    # Verify consistency
    if abs(result_total_reads - total_reads) > total_reads * 0.01:  # allow a 1% tolerance
        print(f"   ⚠️  Warning: the total reads in the result ({result_total_reads:,}) do not match the total input reads ({total_reads:,})")
    
    return result_df


def main():
    parser = argparse.ArgumentParser(
        description='Calculate viral abundance metrics (RPM and RPKM)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # Single sample
  python calculate_abundance_en.py \\
    --bracken results/bracken/sample1_bracken.tsv \\
    --kraken results/kraken2/sample1.report \\
    --output sample1_abundance.tsv
  
  # Batch processing
  for sample in sample1 sample2 sample3; do
    python calculate_abundance_en.py \\
      --bracken results/bracken/${sample}_bracken.tsv \\
      --kraken results/kraken2/${sample}.report \\
      --output ${sample}_abundance.tsv
  done
        """
    )
    
    parser.add_argument('--bracken', required=True,
                       help='Bracken output file path (*.tsv)')
    parser.add_argument('--kraken', required=True,
                       help='Kraken2 report file path (*.report)')
    parser.add_argument('--output', required=True,
                       help='Output file path')
    
    args = parser.parse_args()
    
    print("="*60)
    print("🧬 Viral Abundance Calculator (RPM)")
    print("="*60)
    
    # Read input files
    bracken_df = parse_bracken_output(args.bracken)
    kraken_df = parse_kraken_report(args.kraken)
    
    # Compute the relative abundance (RPM)
    result_df = calculate_rpm(bracken_df, kraken_df)
    
    # Save results
    print(f"\n💾 Saving results to: {args.output}")
    result_df.to_csv(args.output, sep='\t', index=False)
    
    # Print a summary of the statistics
    print("\n" + "="*60)
    print("📈 Statistical Summary")
    print("="*60)
    print(f"Detected species count: {len(result_df)}")
    print(f"Total reads: {result_df['Assigned_Reads'].sum():,}")
    print(f"\nTop 10 most abundant species (by RPM):")
    print("-"*60)
    
    top10 = result_df.head(10)
    for idx, row in top10.iterrows():
        print(f"{row['Species'][:50]:50s} RPM: {row['RPM']:10.2f}  Reads: {row['Assigned_Reads']:>12,}")
    
    print("\n✅ Complete!")


if __name__ == '__main__':
    main()

