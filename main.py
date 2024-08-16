import subprocess
import sys

def run_script(script_name, folder_name):
    try:
        # Run the script with the folder name as an argument
        result = subprocess.run(['python3', script_name, folder_name], check=True)
        print(f"{script_name} executed successfully with folder '{folder_name}'.")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while executing {script_name}: {e}")
        sys.exit(1)

if __name__ == "__main__":
    # Check if a folder name is provided
    if len(sys.argv) < 2:
        print("Please provide the folder name as an argument.")
        sys.exit(1)

    # Get the folder name from the command-line argument
    folder_name = sys.argv[1]

    # Define the scripts to be executed in sequence
    scripts = ["before_trans_mvc.py", "transf_mvc.py", "after_trans_mvc.py","marching_tr.py"]

    # Execute each script with the same folder
    for script in scripts:
        run_script(script, folder_name)

    print("All scripts executed successfully!")
