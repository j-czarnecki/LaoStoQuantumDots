#!/bin/bash

BASE_DIR="$1"

if [ -z "$BASE_DIR" ]; then
    echo "Usage: $0 <base_directory>"
    exit 1
fi

for subdir in "$BASE_DIR"/*/; do
    # skip if no directories found
    [ -d "$subdir" ] || continue

    sbatch <<EOF
#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH -p cpu
#SBATCH --job-name=QD_POST
#SBATCH --output="${subdir}/output.out"
#SBATCH --error="${subdir}/error.err"

srun /home/czarnecki/LaoStoQuantumDots/bin/post_qd.x "$subdir"
EOF

done