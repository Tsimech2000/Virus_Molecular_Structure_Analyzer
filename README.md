# Virus Molecular Structure Analyzer

This Streamlit app allows users to analyze virus molecular structures using data from multiple sources, including RCSB PDB, AlphaFold, PDBe, and SwissModel. Users can upload a PDB file or search by virus name to fetch and visualize molecular structures.

## Features
- **🔍 Search for virus structures** via UniProt and other databases.
- **📂 Upload a PDB file** for analysis.
- **🧬 Extract amino acid sequences** from PDB files.
- **🔬 Render 3D molecular structures** using py3Dmol.
- **📥 Download PDB files** for further analysis.

## Installation
1. Clone the repository:
   ```sh
   git clone https://github.com/yourusername/virus-analyzer.git
   cd virus-analyzer
   ```
2. Install dependencies:
   ```sh
   pip install -r requirements.txt
   ```
3. Run the app:
   ```sh
   streamlit run app.py
   ```

## Deployment on Streamlit Cloud
1. Upload your project to a GitHub repository.
2. Go to [Streamlit Cloud](https://share.streamlit.io/) and connect your GitHub repository.
3. Select `app.py` as the main script.
4. Deploy and share your app!

## License
This project is open-source and free to use. Contributions are welcome!
