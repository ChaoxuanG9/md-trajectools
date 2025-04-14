# md-trajectools

A lightweight collection of Python scripts for post-processing molecular dynamics (MD) simulation trajectories.  
This is **not a Python package**. Instead, this repository provides **ready-to-use analysis scripts** that can be:

- Run directly from the command line
- Copied into your own Jupyter notebooks or workflows
- Modified freely to suit your system and analysis needs

---

## 🔍 What's Included

- `msd.py`: Mean Square Displacement
- `rdf.py`: Radial Distribution Function
- `residence_time.py`: Residence time analysis
- `utils.py`: Shared utility functions
- `examples/`: Usage demos and templates

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

## 📜 License
This project is licensed under the MIT License (see `LICENSE` file). You are free to use, modify, and share this code with attribution.

## 📣 Citation
If you use this code in your research, presentations, or derivative works, please cite this repo.

## 🤝 Contributions
Feedback, bug reports, or improvements are welcome. Feel free to open an issue or submit a pull request!
