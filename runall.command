TITLE='\033[1;95m'  # Formats text pink and bold
NC='\033[0m'        # Removes text formatting

# Ask for confirmation before starting process
read -p "Delete previous results and run full analysis? (y/n) " choice
case $choice in
    [yY]* ) echo "Proceeding..." ;;
    [nN]* ) echo "Cancelled"
    exit ;;
    *) echo "Input not recognized"
    exit ;;
esac


echo "\n${TITLE}~~~~~~~~~~~ ELFIN ANALYSIS BEGIN ~~~~~~~~~~~${NC}"

cd $(dirname -- $0) 			       # Move into folder where this script is located (top level of project folder)
source ./code/python_environment/bin/activate  # Activate virtual Python environment

# Clear out old figures and results to start fresh each time. Existing ELFIN data (raw and processed) is deleted and downloaded in elfin_data_download.py
rm -f ./results/figures/*
echo "rm -f ./results/figures/*"

rm -f ./results/*
echo "rm -f ./results/*"

python3.9 ./code/elfin_data_download.py -delete  # Run Python data retrieval script
rm -rf ./code/__pycache__                	 # Delete __pycache__ directory. Doesn't make any functional difference, just keeps the directory clean and readable.
PYTHON_EXIT_CODE=$?                              # Exit code 0 = all good, exit code 1 = terminated with error

if [ $PYTHON_EXIT_CODE -eq 1 ]; then     # If python exited with error, don't run Julia analysis program
    exit ;;
fi

julia ./code/create_emic_list.jl # Creates EIP event list

julia ./code/create_figures.jl # Generates poster figures

echo "${TITLE}~~~~~~~~~~~ ELFIN ANALYSIS END ~~~~~~~~~~~${NC}\n\n"