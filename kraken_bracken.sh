#!/bin/bash

# Set the script to exit immediately on any error
set -e

# Export paths for Kraken and Bracken
export PATH=/programs/Bracken-2.9:/programs/Bracken-2.9/src:/programs/Bracken-2.9/analysis_scripts:$PATH
export PATH=/programs/kraken2.1.3:$PATH

# Function to display usage
usage() {
    echo "Usage: $0 --fastq-dir DIR [--kraken-db DIR] [--kraken-confidence VALUE] [--kraken-min-hit-groups VALUE] [--threads VALUE] [--bracken] [--bracken-read-length VALUE]"
    exit 1
}

# Parse arguments
FASTQ_DIR=""
KRAKEN_DB=""
KRAKEN_CONFIDENCE=0.1
KRAKEN_MIN_HIT_GROUPS=2
THREADS=1
RUN_BRACKEN=false
BRACKEN_READ_LENGTH=75

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --fastq-dir)
            FASTQ_DIR="$2"
            shift 2
            ;;
        --kraken-db)
            KRAKEN_DB="$2"
            shift 2
            ;;
        --kraken-confidence)
            KRAKEN_CONFIDENCE="$2"
            shift 2
            ;;
        --kraken-min-hit-groups)
            KRAKEN_MIN_HIT_GROUPS="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --bracken)
            RUN_BRACKEN=true
            shift
            ;;
        --bracken-read-length)
            BRACKEN_READ_LENGTH="$2"
            shift 2
            ;;
        *)
            usage
            ;;
    esac
done

# Check if fastq directory is provided
if [[ -z "$FASTQ_DIR" ]]; then
    usage
fi

# Ensure the fastq directory exists
if [[ ! -d "$FASTQ_DIR" ]]; then
    echo "Error: Directory $FASTQ_DIR does not exist."
    exit 1
fi

# Check if kraken database directory is provided
if [[ -z "$KRAKEN_DB" ]]; then
    echo "Error: Kraken database directory must be specified with --kraken-db."
    exit 1
fi

# Ensure the kraken database directory exists
if [[ ! -d "$KRAKEN_DB" ]]; then
    echo "Error: Directory $KRAKEN_DB does not exist."
    exit 1
fi

# Create subdirectories for output organization
KRAKEN_OUTPUT_DIR="$FASTQ_DIR/kraken_output"
KRAKEN_REPORTS_DIR="$FASTQ_DIR/kraken_reports"
BRACKEN_CLASSIFICATION_DIR="$FASTQ_DIR/bracken_classification"
mkdir -p "$KRAKEN_OUTPUT_DIR"
mkdir -p "$KRAKEN_REPORTS_DIR"
mkdir -p "$BRACKEN_CLASSIFICATION_DIR"
echo "Created directories: $KRAKEN_OUTPUT_DIR, $KRAKEN_REPORTS_DIR, $BRACKEN_REPORTS_DIR, $BRACKEN_CLASSIFICATION_DIR"

# Process each FASTQ file in the directory
for file in "$FASTQ_DIR"/*.fastq; do
    if [[ -f "$file" ]]; then
        base=$(basename "$file" .fastq)
        echo "Processing file: $file"
        kraken2 \
            --db "$KRAKEN_DB" \
            --confidence "$KRAKEN_CONFIDENCE" \
            --minimum-hit-groups "$KRAKEN_MIN_HIT_GROUPS" \
            --threads "$THREADS" \
            --report "$KRAKEN_REPORTS_DIR/${base}_report.txt" \
            --output "$KRAKEN_OUTPUT_DIR/${base}_output.txt" \
            "$file"
        echo "Finished processing: $file"
    else
        echo "Skipping non-file: $file"
    fi
done

# Combine and format Kraken report files
combine_reports() {
    REPORT_DIR="$KRAKEN_REPORTS_DIR"

    # Define the column headers with "%_reads"
    headers="SampleID\t%_reads\tTotal_Reads\tClassified_Reads\tRankCode\tNCBITaxID\tTaxa"

    # Create a temporary file to store the combined output
    output_file=$(mktemp)
    echo -e "$headers" > "$output_file"

    # Process each _report.txt file in the directory
    for file in "$REPORT_DIR"/*_report.txt; do
        if [ -f "$file" ]; then
            filename=$(basename "$file")

            # Extract the sample ID as everything before "_report.txt"
            sample_id="${filename%%_report.txt}"

            echo "Formatting file: $filename"

            # Create a temporary file for the formatted content
            tmp_file=$(mktemp)

            # Add headers to the temporary file
            echo -e "$headers" > "$tmp_file"

            # Format the file and ensure the first row has the SampleID
            awk -v sample_id="$sample_id" 'BEGIN {OFS="\t"} {
                if (NR == 1) {
                    # Add SampleID to the first row and shift other fields one to the right
                    print sample_id, $0;
                } else {
                    # Process subsequent rows as usual
                    print sample_id, $0;
                }
            }' "$file" >> "$tmp_file"

            # Append the cleaned file's content to the combined output
            tail -n +2 "$tmp_file" >> "$output_file"
        fi
    done

    # Save the combined content to the final output file
    final_output_file="$REPORT_DIR/combined_report.txt"
    mv "$output_file" "$final_output_file"

    echo "Combined report saved to: $final_output_file"

    # Split the combined report into taxonomic ranks
    split_taxa_ranks "$final_output_file" "$REPORT_DIR"
}

# Split taxonomic ranks
split_taxa_ranks() {
    combined_file="$1"
    output_dir="$2"

    # Validate the combined file
    if [ ! -f "$combined_file" ]; then
        echo "Error: The file '$combined_file' does not exist."
        exit 1
    fi

    # Define taxonomic ranks and their file names
    declare -A rank_files=(
        ["U"]="combined_report_unclassified.txt"
        ["R"]="combined_report_root.txt"
        ["D"]="combined_report_domain.txt"
        ["K"]="combined_report_kingdom.txt"
        ["P"]="combined_report_phylum.txt"
        ["C"]="combined_report_class.txt"
        ["O"]="combined_report_order.txt"
        ["F"]="combined_report_family.txt"
        ["G"]="combined_report_genus.txt"
        ["S"]="combined_report_species.txt"
    )

    # Initialize empty files for each rank
    for rank in "${!rank_files[@]}"; do
        echo -e "SampleID\t%_reads\tTotal_Reads\tClassified_Reads\tRankCode\tNCBITaxID\tTaxa" > "$output_dir/${rank_files[$rank]}"
    done

    # Process the combined file and separate by taxonomic rank
    tail -n +2 "$combined_file" | while IFS=$'\t' read -r SampleID reads Total_Reads Classified_Reads RankCode NCBITaxID Taxa; do
        # Extract the base rank (e.g., "G1" -> "G", "S2" -> "S")
        base_rank="${RankCode:0:1}"

        # Check if the rank matches any predefined ranks
        if [[ -v rank_files[$base_rank] ]]; then
            # Append the row to the corresponding file
            echo -e "$SampleID\t$reads\t$Total_Reads\t$Classified_Reads\t$RankCode\t$NCBITaxID\t$Taxa" >> "$output_dir/${rank_files[$base_rank]}"
        fi
    done

    echo "Taxonomic ranks have been separated into files in the '$output_dir' directory."
}

# Run Bracken
run_bracken() {
    TEMP_KRAKEN_REPORTS_DIR="$KRAKEN_REPORTS_DIR/temp_kraken_reports"
    mkdir -p "$TEMP_KRAKEN_REPORTS_DIR"

    echo "Moving _report.txt files to a temporary directory..."
    for file in "$KRAKEN_REPORTS_DIR"/*_report.txt; do
        if [[ -f "$file" ]]; then
            mv "$file" "$TEMP_KRAKEN_REPORTS_DIR/"
        fi
    done

    for report in "$TEMP_KRAKEN_REPORTS_DIR"/*_report.txt; do
        if [[ -f "$report" ]]; then
            base=$(basename "$report" _report.txt)
            echo "Running Bracken for file: $report..."

            bracken -d "$KRAKEN_DB" \
                -i "$report" \
                -o "$BRACKEN_CLASSIFICATION_DIR/${base}_bracken_species.txt" \
                -l S \
                -r "$BRACKEN_READ_LENGTH"

            bracken -d "$KRAKEN_DB" \
                -i "$report" \
                -o "$BRACKEN_CLASSIFICATION_DIR/${base}_bracken_S1.txt" \
                -l S1 \
                -r "$BRACKEN_READ_LENGTH"

            bracken -d "$KRAKEN_DB" \
                -i "$report" \
                -o "$BRACKEN_CLASSIFICATION_DIR/${base}_bracken_S2.txt" \
                -l S2 \
                -r "$BRACKEN_READ_LENGTH"

            echo "Bracken files generated for $report"
        fi
    done

    echo "Restoring _report.txt files to $KRAKEN_REPORTS_DIR..."
    for file in "$TEMP_KRAKEN_REPORTS_DIR"/*_report.txt; do
        if [[ -f "$file" ]]; then
            mv "$file" "$KRAKEN_REPORTS_DIR/"
        fi
    done

    echo "Cleaning up temporary directory..."
    rm -rf "$TEMP_KRAKEN_REPORTS_DIR"
}

# Combine Bracken Outputs
combine_bracken_outputs() {
    INPUT_DIR="$1"
    OUTPUT_FILE="$INPUT_DIR/bracken_combined.txt"

    if ! ls "$INPUT_DIR"/*_bracken_*.txt 1> /dev/null 2>&1; then
        echo "Error: No _bracken_*.txt files found in $INPUT_DIR."
        exit 1
    fi

    if [ -f "$OUTPUT_FILE" ]; then
        rm "$OUTPUT_FILE"
    fi

    HEADER_WRITTEN=false

    for FILE in "$INPUT_DIR"/*_bracken_*.txt; do
        SAMPLE_ID=$(basename "$FILE" | sed -E 's/_bracken.*//')

        while IFS= read -r LINE; do
            if [[ "$LINE" == "name"* && "$HEADER_WRITTEN" == "false" ]]; then
                echo -e "${LINE}\tSampleID" >> "$OUTPUT_FILE"
                HEADER_WRITTEN=true
            elif [[ "$LINE" != "name"* ]]; then
                echo -e "${LINE}\t${SAMPLE_ID}" >> "$OUTPUT_FILE"
            fi
        done < "$FILE"
    done

    echo "Combined file created: $OUTPUT_FILE"
}

# Run Bracken
run_bracken

# Combine Bracken Outputs
combine_bracken_outputs "$BRACKEN_CLASSIFICATION_DIR"

echo "Pipeline completed successfully. Outputs are in $KRAKEN_OUTPUT_DIR, $KRAKEN_REPORTS_DIR, and $BRACKEN_CLASSIFICATION_DIR."
