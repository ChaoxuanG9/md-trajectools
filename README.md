# md-trajectools

A lightweight collection of Python scripts for post-processing molecular dynamics (MD) simulation trajectories.  
This is **not a Python package**. Instead, this repository provides ready-to-use and reusable analysis scripts that can be:

- Run directly from the command line
- Copied into your own Jupyter notebooks or workflows
- Modified freely to suit your system and analysis needs

---

## 🔍 What's Included

The functionalities included are:

- **Residence Time**: Measure how long a molecule or ion stays within a defined region.
- **Mean Squared Displacement (MSD)**: Analyze particle diffusion, either the "absolute" MSD from time zero or the time-averaged MSD, with the option of correcting for the movement of the system's center of mass.
- **Plotting the Electric Double Layer (EDL) Solution Environment**: Visualize the distribution of species in the EDL based on the ion solvation environment.
- **Charge density, electric field, and electrostatic potential**: Calculate all-atom charge density, electric field intensity, and electrostatic potential drop in slab-geometry MD simulations.
- **Computing the Radial Distribution Function (RDF)**: Calculate the probability of finding a particle at a certain distance from a reference particle.

More to be updated...

---

## 💡 Use Philosophy
This project follows a simple principle: "Transparent, editable, and reusable". You are encouraged to:
- Run these scripts as-is
- Copy-paste functions into your own workflow
- Modify for your specific MD system or output format

---

## ⚡ Quick Start

```bash
# Clone the repo
git clone https://github.com/yourusername/md-trajectools.git
cd md-trajectools

# Install required dependencies
pip install -r requirements.txt

# Run an analysis script directly
python msd.py --traj trajectory.xyz
```
---

## 📂 Repository Structure

```
md-trajectools/
├── README.md                  # Clear project intro, instructions
├── LICENSE                    # Open-source clarity
├── residence_time/            # Logical grouping by function
│   ├── residence_time.py
│   └── config_template.yaml 
├── MSD/                       # MSD-specific tools
│   └── msd.py
├── EDL_structure/             # EDL plotting scripts
│   └── EDL_plot.py
├── Electric_potential/        # Potential profile computation
│   └── electric_potential.py
├── RDF/                       # RDF calculators
│   └── rdf.py
├── docs/                      # For future documentation or Sphinx
├── examples/                  # Ready-to-run usage demos
│   ├── example_trajectory.xyz
│   └── example_usage.sh


```

- `residence_time.py`: Script for calculating residence time.
- `MSD.py`: Script for calculating mean squared displacement with center of mass drift correction.
- `EDL_plot.py`: Script for plotting the electric double layer solution environment.
- `electric_potential.py`: Script for calculating all-atom charge density, electric field intensity, and electrostatic potential of slab-geometry MD simulations.
- `RDF.py`: Script for computing the radial distribution function.
- `examples/`: Directory containing example trajectory files along with example usage script.

---
## 📜 License
This project is licensed under the MIT License (see `LICENSE` file). You are free to use, modify, and share this code with attribution.

## 📣 Citation
If you use this code in your research, presentations, or derivative works, please cite this repo.

## 🤝 Contributions
Feedback, bug reports, or improvements are welcome. Feel free to open an issue or submit a pull request!
