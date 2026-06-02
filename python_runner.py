import subprocess
import sys

command = ["nextflow", "run", "main.nf"]

# Append Docker runtime arguments
command.extend(sys.argv[1:])

print("Running command:")
print(" ".join(command))

process = subprocess.Popen(
    command,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
    text=True
)

for line in process.stdout:
    print(line, end="")

process.wait()

print(f"\nProcess finished with code {process.returncode}")
