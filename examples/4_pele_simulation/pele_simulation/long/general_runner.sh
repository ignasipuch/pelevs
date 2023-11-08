for d in *; do if [ "$d" != "general_runner.sh" ]; then cd "$d"; sbatch run_plat; cd ..; fi; done
