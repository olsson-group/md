for i in $(ls $1); do
    echo "Simulating $1/$i"
    python slurm/submit.py \"python scripts/simulate.py $1/$i/traj.pdb\"
done
