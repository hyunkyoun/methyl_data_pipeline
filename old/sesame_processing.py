import subprocess
import numpy as np
import json

def call_script(file_path):
    result = subprocess.run(['Rscript', 'sesame_processing.R', file_path],
                            capture_output=True,
                            text=True,
                            check=True
    )

    matrix = np.array(json.loads(result.stdout))
    return matrix

if __name__ == "__main__":
    file_path = './data/idat_files'
    matrix = call_script(file_path)
    print(matrix)