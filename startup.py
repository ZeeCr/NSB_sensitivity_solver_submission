import sys
import csv
import os
import subprocess

id_idx_file = './data/id_idx.dat'
with open(id_idx_file, 'r') as file:
    id_idx = int(file.read())

try:
    subprocess.run([f"screen -dmS {id_idx} -L"], shell = True)
    subprocess.run([f"screen -r {id_idx}"], shell = True)
except KeyboardInterrupt:
    subprocess.run(["echo", "Keyboard interrupt detected. Exiting..."])
    sys.exit(0)