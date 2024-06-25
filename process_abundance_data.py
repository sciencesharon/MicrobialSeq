import pandas as pd
import os
from pathlib import Path
import argparse

def merged_abundance_processing(merged, output_dir):
    # Step 1: Make a temporary tax table
    def merged_abundance_to_tax(merged):
        tax = merged.iloc[:, [0]].copy()
        tax = tax.iloc[2:, :].reset_index(drop=True)
        tax.index = ['OTU' + str(i) for i in range(len(tax))]
        tax.columns = ['V1']
        tax = tax['V1'].str.split('|', expand=True)
        tax.columns = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"]
        
        for col in tax.columns:
            tax[col] = tax[col].str.replace(f'{col[0].lower()}__', '', regex=True)
            tax[col] = tax[col].str.replace('_', ' ', regex=True)
        
        return tax

    # Step 2: Make a temporary OTU table
    def merged_abundance_to_otu(merged):
        merged.columns = merged.iloc[0]
        merged = merged.drop([0, 1]).reset_index(drop=True)
        merged = merged[~merged['clade_name'].str.contains('Archaea|Viruses|Eukaryota')]
        otu = merged.drop(columns=['clade_name']).apply(pd.to_numeric)
        otu.index = ['OTU' + str(i) for i in range(len(otu))]
        return otu

    # Step 3: Subset the tax table by !is.na for each column and save as a new dataframe
    def subset_by_na(df):
        df_list = {col_name: df[df[col_name].notna()][[col_name]] for col_name in df.columns}
        return df_list

    # Step 4: Subset OTU table by rownames of each new tax table
    def subset_by_rownames(df_list, dataframe2):
        df_subset_list = {df_name: dataframe2.loc[dataframe2.index.isin(df_list[df_name].index)] for df_name in df_list}
        return df_subset_list

    # Step 5: Now we combine the OTU and tax tables
    def subset_and_cbind(df_list, dataframe2):
        df_combined_list = {df_name: pd.concat([df_list[df_name], dataframe2.loc[df_list[df_name].index]], axis=1) for df_name in df_list}
        return df_combined_list

    # Step 6: Combine rows which are identical taxa
    def summarize_and_save_to_csv(df_list, output_dir):
        summarized_list = {}
        for df_name, df in df_list.items():
            df = df[~df.iloc[:, 0].str.contains('Unclassified|GB|unclassified', case=False)]
            unique_values = df.iloc[:, 0].unique()
            result_df = pd.DataFrame()
            for value in unique_values:
                subset_df = df[df.iloc[:, 0] == value]
                summed_row = subset_df.iloc[:, 1:].sum(axis=0)
                new_row = pd.DataFrame([summed_row], index=[value])
                result_df = pd.concat([result_df, new_row], axis=0)
            result_df.index.name = df.columns[0]
            new_df_name = f"{df_name}_ReadCount"
            result_df.to_csv(Path(output_dir) / f"{new_df_name}.csv")
            summarized_list[new_df_name] = result_df
        return summarized_list

    # Step 7: Convert to Relative Abundance
    def convert_to_percentages_list(df_list, output_dir):
        df_percent_list = {}
        for df_name, df in df_list.items():
            df_percent = df.apply(lambda x: (x / x.sum()) * 100)
            new_df_name = f"{df_name.split('_')[0]}_RelAbun"
            df_percent.to_csv(Path(output_dir) / f"{new_df_name}.csv")
            df_percent_list[new_df_name] = df_percent
        return df_percent_list

    tax = merged_abundance_to_tax(merged)
    otu = merged_abundance_to_otu(merged)
    df_list = subset_by_na(tax)
    df_list2 = subset_by_rownames(df_list, otu)
    df_combined_list = subset_and_cbind(df_list, otu)
    summarized_list = summarize_and_save_to_csv(df_combined_list, output_dir)
    percentage_list = convert_to_percentages_list(summarized_list, output_dir)

    return percentage_list

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process merged abundance data and save results as CSV files.")
    parser.add_argument('input_file', type=str, help="Path to the input merged CSV file.")
    parser.add_argument('output_dir', type=str, help="Directory where the output CSV files will be saved.")

    args = parser.parse_args()

    # Load the merged data
    merged = pd.read_csv(args.input_file)

    # Process the merged data
    processed_list = merged_abundance_processing(merged, args.output_dir)

    print("Processing completed. The results are saved in the specified output directory.")
