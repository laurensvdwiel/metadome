#!/bin/bash
set -e

# Default values
COMPOSE_FILE="docker-compose-populate-database.yml"

# These must be set by command line parameters
DATA_DIR=""
FILE_PATH=""

# Display usage information
function show_usage {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -d, --data-dir DIR     Data directory path"
    echo "  -f, --file FILE        Input file name"
    echo "  -c, --compose-file FILE Compose file to use (default: $COMPOSE_FILE)"
    echo "  -h, --help             Show this help message"
    echo
    echo "Example: $0 --data-dir /path/to/data --file my_input.csv.gz"
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -d|--data-dir)
            DATA_DIR="$2"
            shift 2
            ;;
        -f|--file)
            FILE_PATH="$2"
            shift 2
            ;;
        -c|--compose-file)
            COMPOSE_FILE="$2"
            shift 2
            ;;
        -h|--help)
            show_usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            show_usage
            exit 1
            ;;
    esac
done

# Check required parameters
if [ -z "$DATA_DIR" ]; then
    echo "ERROR: Data directory path (-d, --data-dir) is required."
    show_usage
    exit 1
fi

if [ -z "$FILE_PATH" ]; then
    echo "ERROR: Input file name (-f, --file) is required."
    show_usage
    exit 1
fi


# Validate parameters
if [ ! -d "$DATA_DIR" ]; then
    echo "Error: Data directory '$DATA_DIR' does not exist."
    exit 1
fi

# If FILE_PATH is relative, it will be relative to the container's /prebuild_data
# So we don't need to check if it exists in the host's DATA_DIR
echo "Using data directory: $DATA_DIR"
echo "Using input file: $FILE_PATH"
echo "Using compose file: $COMPOSE_FILE"

# Full command in the container
CONTAINER_FILE_PATH="/prebuild_data/$FILE_PATH"
CONTAINER_CMD="python -m metadome.batch_load --csv $CONTAINER_FILE_PATH"

echo "Running: $CONTAINER_CMD"
echo "Starting application... This may take a while."

# Run the command
docker compose -f "$COMPOSE_FILE" run --rm \
  -v "$DATA_DIR:/prebuild_data" \
  app $CONTAINER_CMD

# Exit code of the previous command
EXIT_CODE=$?

# Always clean up, regardless of success or failure
echo "Application finished with exit code: $EXIT_CODE"
echo "Cleaning up environment..."
docker compose -f "$COMPOSE_FILE" down \
  --remove-orphans \
  -v \
  --rmi local

echo "Cleanup complete."

# Return the original exit code from the application
exit $EXIT_CODE
