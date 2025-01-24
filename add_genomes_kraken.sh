#!/bin/bash

# Exit script on any error
set -e

# Function to display usage
usage() {
    echo "Usage: $0 [--builddb] [--db-dir <path>] [--db-threads <num>] [--fna-dir <path>] [--python-script <path>] [--kraken-path <path>] [--bracken-path <path>] [--update-fna] [--append-library] [--append-prelim] [--append-names] [--append-nodes] [--cleanup] [--rebuild-db] [--build-bracken] [--read-length <num>] [--kmer-length <num>] [--strain-nodes-file <path>] [--help]"
    echo "Options:"
    echo "  --builddb               Build a Kraken2 database with bacteria"
    echo "  --db-dir <path>         Specify an existing or new directory for the Kraken2 database (default: ./kraken2db)"
    echo "  --db-threads <num>      Specify the number of threads to use for building the database (default: 4)"
    echo "  --fna-dir <path>        Directory containing .fna files and strain info .txt file (required for --update-fna, --append-library, --append-prelim, --append-names, and --append-nodes)"
    echo "  --python-script <path>  Specify the Python script to use for updating .fna headers (default: /workdir/sct88/shells/updatefna.py)"
    echo "  --kraken-path <path>    Specify the path to the Kraken2 installation (default: /programs/kraken2.1.3)"
    echo "  --bracken-path <path>   Specify the path to the Bracken installation (default: /programs/Bracken-2.9)"
    echo "  --update-fna            Run the Python script to update .fna headers"
    echo "  --append-library        Append cleaned .fna files to the Kraken2 library"
    echo "  --append-prelim         Append cleaned .fna files to prelim_map.txt"
    echo "  --append-names          Append cleaned .fna files to names.dmp"
    echo "  --append-nodes          Append entries to nodes.dmp using parent and strain taxids"
    echo "  --cleanup               Remove .k2d, .map, and .txt files from the top-level Kraken2 database directory"
    echo "  --rebuild-db            Rebuild the Kraken2 database after cleanup"
    echo "  --build-bracken         Build Bracken files after rebuilding the Kraken2 database"
    echo "  --read-length <num>     Specify the read length for Bracken build (required for --build-bracken)"
    echo "  --kmer-length <num>     Specify the k-mer length for Bracken build (required for --build-bracken)"
    echo "  --strain-nodes-file <path>  Specify the path to the parent_and_strain_ids.txt file (optional if --fna-dir is provided)"
    echo "  --help                  Display this help message"
    exit 1
}

# Default parameters
DB_DIR="./kraken2db"
DB_THREADS=1
BUILD_DB=false
UPDATE_FNA=false
APPEND_LIBRARY=false
APPEND_PRELIM=false
APPEND_NAMES=false
APPEND_NODES=false
CLEANUP=false
REBUILD_DB=false
BUILD_BRACKEN=false
READ_LENGTH=""
KMER_LENGTH=""
FNA_DIR=""
STRAIN_NODES_FILE=""
PYTHON_SCRIPT="/workdir/sct88/shells/updatefna.py"
KRAKEN_PATH="/programs/kraken2.1.3"
BRACKEN_PATH="/programs/Bracken-2.9"

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --builddb) BUILD_DB=true ;;
        --db-dir) DB_DIR="$2"; shift ;;
        --db-threads) DB_THREADS="$2"; shift ;;
        --fna-dir) FNA_DIR="$2"; shift ;;
        --python-script) PYTHON_SCRIPT="$2"; shift ;;
        --kraken-path) KRAKEN_PATH="$2"; shift ;;
        --bracken-path) BRACKEN_PATH="$2"; shift ;;
        --update-fna) UPDATE_FNA=true ;;
        --append-library) APPEND_LIBRARY=true ;;
        --append-prelim) APPEND_PRELIM=true ;;
        --append-names) APPEND_NAMES=true ;;
        --append-nodes) APPEND_NODES=true ;;
        --cleanup) CLEANUP=true ;;
        --rebuild-db) REBUILD_DB=true ;;
        --build-bracken) BUILD_BRACKEN=true ;;
        --read-length) READ_LENGTH="$2"; shift ;;
        --kmer-length) KMER_LENGTH="$2"; shift ;;
        --strain-nodes-file) STRAIN_NODES_FILE="$2"; shift ;;
        --help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
    shift
done

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --builddb) BUILD_DB=true ;;
        --db-dir) DB_DIR="$2"; shift ;;
        --db-threads) DB_THREADS="$2"; shift ;;
        --fna-dir) FNA_DIR="$2"; shift ;;
        --update-fna) UPDATE_FNA=true ;;
        --append-library) APPEND_LIBRARY=true ;;
        --append-prelim) APPEND_PRELIM=true ;;
        --append-names) APPEND_NAMES=true ;;
        --append-nodes) APPEND_NODES=true ;;
        --cleanup) CLEANUP=true ;;
        --rebuild-db) REBUILD_DB=true ;;
        --build-bracken) BUILD_BRACKEN=true ;;
        --read-length) READ_LENGTH="$2"; shift ;;
        --kmer-length) KMER_LENGTH="$2"; shift ;;
        --strain-nodes-file) STRAIN_NODES_FILE="$2"; shift ;;
        --help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
    shift
done

# Update PATH for Kraken2
export PATH="$KRAKEN_PATH:$PATH"
echo "Kraken2 path set to: $KRAKEN_PATH"

# Ensure Kraken2 is available
if ! command -v kraken2-build &> /dev/null; then
    echo "Error: kraken2-build command not found in $KRAKEN_PATH. Ensure Kraken2 is installed."
    exit 1
fi

# Update PATH for Bracken
export PATH="$BRACKEN_PATH:$BRACKEN_PATH/src:$BRACKEN_PATH/analysis_scripts:$PATH"
echo "Bracken path set to: $BRACKEN_PATH"

# Ensure Bracken is available if building Bracken files
if [[ "$BUILD_BRACKEN" == true ]]; then
    if ! command -v bracken-build &> /dev/null; then
        echo "Error: bracken-build command not found in $BRACKEN_PATH. Ensure Bracken is installed."
        exit 1
    fi
    if [[ -z "$READ_LENGTH" ]]; then
        echo "Error: --read-length must be specified for --build-bracken."
        exit 1
    fi
    if [[ -z "$KMER_LENGTH" ]]; then
        echo "Error: --kmer-length must be specified for --build-bracken."
        exit 1
    fi
fi

# Step 1: Build the Kraken2 database (if requested)
if [[ "$BUILD_DB" == true ]]; then
    echo "Building Kraken2 database..."
    mkdir -p "$DB_DIR"  # Ensure the database directory exists
    echo "Using database directory: $DB_DIR"
    kraken2-build --download-library bacteria --db "$DB_DIR" --threads "$DB_THREADS"
    kraken2-build --build --db "$DB_DIR" --threads "$DB_THREADS"
    echo "Kraken2 database built successfully in $DB_DIR."
else
    echo "Using existing database directory: $DB_DIR"
    if [[ ! -d "$DB_DIR" ]]; then
        echo "Error: Specified database directory $DB_DIR does not exist."
        exit 1
    fi
fi

# Step 2: Update FNA headers (if requested)
if [[ "$UPDATE_FNA" == true ]]; then
    if [[ -z "$FNA_DIR" ]]; then
        echo "Error: --fna-dir must be specified to update .fna headers."
        exit 1
    fi
    if [[ ! -f "$PYTHON_SCRIPT" ]]; then
        echo "Error: Specified Python script not found at $PYTHON_SCRIPT."
        exit 1
    fi
fi
    STRAIN_INFO_FILE="$FNA_DIR/straininfo.txt"
    if [[ ! -f "$STRAIN_INFO_FILE" ]]; then
        echo "Error: Strain information file not found in $FNA_DIR. Expected file: straininfo.txt."
        exit 1
    fi

    echo "Running Python script to update .fna headers..."
    python3 "$PYTHON_SCRIPT" "$STRAIN_INFO_FILE" "$FNA_DIR"
    echo "Headers updated. Cleaned files are in $FNA_DIR/cleaned."
fi

# Step 3: Append cleaned .fna files to the Kraken2 library (if requested)
if [[ "$APPEND_LIBRARY" == true ]]; then
    LIBRARY_FNA="$DB_DIR/library/bacteria/library.fna"
    if [[ ! -d "$DB_DIR/library/bacteria" ]]; then
        echo "Error: Bacteria library directory not found in $DB_DIR. Ensure the database is built or correctly specified."
        exit 1
    fi

    echo "Appending cleaned .fna files to the Kraken2 library..."
    mkdir -p "$(dirname "$LIBRARY_FNA")"
    cat "$FNA_DIR/cleaned/"*.fna >> "$LIBRARY_FNA"
    echo "Cleaned .fna files appended to $LIBRARY_FNA."
fi

# Step 4: Append to prelim_map.txt (if requested)
if [[ "$APPEND_PRELIM" == true ]]; then
    PRELIM_MAP="$DB_DIR/library/bacteria/prelim_map.txt"
    if [[ ! -d "$DB_DIR/library/bacteria" ]]; then
        echo "Error: Bacteria library directory not found in $DB_DIR. Ensure the database is built or correctly specified."
        exit 1
    fi

    echo "Appending to prelim_map.txt..."
    awk '/^>/ {
        gsub("^>", "");
        split($1, id, "|");
        entry = "TAXID\tkraken:taxid|" id[3] "|" id[1] "\t" id[3];
        if (!(entry in seen)) {
            print entry;
            seen[entry] = 1;
        }
    }' "$FNA_DIR/cleaned/"*.fna >> "$PRELIM_MAP"
    echo "Entries from cleaned .fna files appended to $PRELIM_MAP."
fi

# Step 5: Append to names.dmp (if requested)
if [[ "$APPEND_NAMES" == true ]]; then
    NAMES_DMP="$DB_DIR/taxonomy/names.dmp"
    if [[ ! -d "$DB_DIR/taxonomy" ]]; then
        echo "Error: Taxonomy directory not found in $DB_DIR. Ensure the database is built or correctly specified."
        exit 1
    fi

    echo "Appending to names.dmp..."
    awk '/^>/ {
        gsub("^>", "");
        split($1, id, "|");
        print id[3] "\t|\t" id[1] "\t|\t\t|\tscientific name\t|";
    }' "$FNA_DIR/cleaned/"*.fna | sort -u >> "$NAMES_DMP"
    echo "Entries from cleaned .fna files appended to $NAMES_DMP."
fi

# Step 6: Append to nodes.dmp (if requested)
if [[ "$APPEND_NODES" == true ]]; then
    if [[ -z "$STRAIN_NODES_FILE" ]]; then
        STRAIN_NODES_FILE="$FNA_DIR/parent_and_strain_ids.txt"
    fi

    if [[ ! -f "$STRAIN_NODES_FILE" ]]; then
        echo "Error: Parent and strain IDs file not found. Expected file: $STRAIN_NODES_FILE."
        exit 1
    fi

    NODES_DMP="$DB_DIR/taxonomy/nodes.dmp"
    if [[ ! -d "$DB_DIR/taxonomy" ]]; then
        echo "Error: Taxonomy directory not found in $DB_DIR. Ensure the database is built or correctly specified."
        exit 1
    fi

    echo "Appending to nodes.dmp..."
    awk '
    BEGIN {
        format = "%d\t|\t%d\t|\tstrain\t|\t\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0|\t1\t|\t0\t|\t0\t|\t\t|\n";
    }
    {
        parent_taxid = $1;  # First column: parent taxonomic ID
        strain_taxid = $2;  # Second column: strain taxonomic ID
        printf(format, strain_taxid, parent_taxid);
    }' "$STRAIN_NODES_FILE" >> "$NODES_DMP"
    echo "Entries from $STRAIN_NODES_FILE appended to $NODES_DMP."
fi

# Step 7: Cleanup (if requested)
if [[ "$CLEANUP" == true ]]; then
    echo "Cleaning up files in the top-level of $DB_DIR..."
    find "$DB_DIR" -maxdepth 1 -type f \( -name "*.k2d" -o -name "*.map" -o -name "*.txt" \) -exec rm -v {} +
    echo "Cleanup completed."
fi

# Step 8: Rebuild Kraken2 database (if requested)
if [[ "$REBUILD_DB" == true ]]; then
    echo "Rebuilding the Kraken2 database in $DB_DIR..."
    kraken2-build --build --db "$DB_DIR" --threads "$DB_THREADS"
    echo "Rebuild of Kraken2 database completed."
fi

# Step 9: Build Bracken files (if requested)
if [[ "$BUILD_BRACKEN" == true ]]; then
    echo "Building Bracken files..."
    bracken-build -d "$DB_DIR" -l "$READ_LENGTH" -k "$KMER_LENGTH" -t "$DB_THREADS"
    echo "Bracken files built successfully in $DB_DIR."
fi

echo "Script completed successfully!"
