import streamlit as st
import os
import tempfile
import subprocess
import py3Dmol
import requests
import time
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from Bio.PDB.Polypeptide import is_aa

# --- Set PyMOL Path (Update this manually) ---
PYMOL_PATH = r"C:\Users\juano\.conda\envs\pymol_env\Scripts\pymol.exe"

# --- Streamlit App Title ---
st.title("🦠 Virus Molecular Structure Analyzer (Multi-Source: RCSB PDB, AlphaFold, PDBe, SwissModel, NCBI MMDB, JGI IMG)")

# --- Sidebar Section ---
st.sidebar.header("📂 Search, Fetch, or Upload a Virus PDB File")
uploaded_file = st.sidebar.file_uploader("📂 Upload a PDB file", type=["pdb"])
virus_name = st.sidebar.text_input("🔍 Enter a virus name (e.g., 'SARS-CoV-2', 'Influenza', 'HIV')")

# --- Function to Get Protein ID from UniProt ---
def get_uniprot_id(virus_name):
    """Fetches the correct protein ID from UniProt based on the virus name."""
    st.sidebar.info("🔍 Searching UniProt for related proteins...")

    url = f"https://rest.uniprot.org/uniprotkb/search?query={virus_name}&format=json"
    response = requests.get(url)

    if response.status_code == 200:
        data = response.json()
        if "results" in data and len(data["results"]) > 0:
            protein_name = data["results"][0]["proteinDescription"]["recommendedName"]["fullName"]["value"]
            uniprot_id = data["results"][0]["primaryAccession"]
            st.sidebar.success(f"✅ Found protein: {protein_name} (UniProt ID: {uniprot_id})")
            return uniprot_id

    st.sidebar.warning("⚠️ No protein match found in UniProt.")
    return None

# --- Function to Search for Structures in Multiple Databases ---
def search_structure(uniprot_id):
    """Searches for protein-related structures in multiple databases."""
    if not uniprot_id:
        return None

    st.sidebar.info("🔍 Searching Structure Databases...")

    sources = {
        "RCSB PDB": f"https://files.rcsb.org/download/{uniprot_id}.pdb",
        "PDBe": f"https://www.ebi.ac.uk/pdbe/entry-files/download/{uniprot_id}.pdb",
        "SwissModel": f"https://swissmodel.expasy.org/repository/uniprot/{uniprot_id}.pdb"
    }

    for source, url in sources.items():
        response = requests.get(url)
        if response.status_code == 200:
            st.sidebar.success(f"✅ Found structure in {source}!")
            return url  # Return valid PDB URL

    st.sidebar.warning("⚠️ No structure found in available databases.")
    return None

# --- Function to Fetch and Save PDB File ---
def fetch_pdb_file(source_url):
    """Downloads a PDB file and saves it locally for analysis."""
    if not source_url:
        return None

    st.sidebar.info("📡 Downloading PDB file...")

    response = requests.get(source_url)
    if response.status_code == 200:
        temp_pdb_path = os.path.join(tempfile.gettempdir(), "downloaded_structure.pdb")
        with open(temp_pdb_path, "wb") as f:
            f.write(response.content)
        st.sidebar.success(f"✅ Successfully downloaded PDB file.")
        return temp_pdb_path

    st.sidebar.warning("❌ Structure download failed.")
    return None

# --- Function to Extract Amino Acid Sequence from Downloaded PDB ---
def extract_sequence_from_pdb(pdb_path):
    """Extracts amino acid sequences from the downloaded PDB file."""
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("molecule", pdb_path)
        sequences = []

        for model in structure:
            for chain in model:
                seq = ""
                for residue in chain:
                    if is_aa(residue, standard=True):  # Ensuring it's a valid amino acid
                        try:
                            seq += seq1(residue.get_resname())  # Converts 3-letter code to 1-letter code
                        except KeyError:
                            seq += "X"  # Placeholder for unknown residues
                if seq:
                    sequences.append(f"Chain {chain.id}: {seq}")

        return "\n".join(sequences) if sequences else "No valid sequence found in PDB."
    
    except Exception as e:
        return f"Error extracting sequence: {e}"

# --- Function to Extract Structural Information from Downloaded PDB ---
def extract_chemical_info(pdb_path):
    """Extracts molecular properties from the downloaded PDB file."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("virus_structure", pdb_path)

    atom_count = 0
    residue_count = 0
    hetatom_count = 0
    chains = set()

    for model in structure:
        for chain in model:
            chains.add(chain.id)
            for residue in chain:
                residue_count += 1
                for atom in residue:
                    atom_count += 1
                    if atom.id[0] == "H":
                        hetatom_count += 1

    mol_weight = atom_count * 12  # Approximate assuming carbon-based backbone

    return {
        "Number of Chains": len(chains),
        "Total Atoms": atom_count,
        "Total Residues": residue_count,
        "Total Hetatoms": hetatom_count,
        "Estimated Molecular Weight (Da)": mol_weight,
    }

# --- Handling PDB File ---
pdb_file = None
sequence = None
chemical_info = None

if uploaded_file:
    pdb_file = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb").name
    with open(pdb_file, "wb") as f:
        f.write(uploaded_file.read())
    sequence = extract_sequence_from_pdb(pdb_file)
    chemical_info = extract_chemical_info(pdb_file)

elif virus_name:
    uniprot_id = get_uniprot_id(virus_name)
    if uniprot_id:
        source_url = search_structure(uniprot_id)
        if source_url:
            pdb_file = fetch_pdb_file(source_url)
            if pdb_file:
                sequence = extract_sequence_from_pdb(pdb_file)
                chemical_info = extract_chemical_info(pdb_file)

# --- Provide Downloadable PDB File ---
if pdb_file and os.path.exists(pdb_file):
    with open(pdb_file, "rb") as f:
        st.sidebar.download_button(
            label="📥 Download PDB File",
            data=f,
            file_name="virus_structure.pdb",
            mime="chemical/x-pdb"
        )

# --- Display Chemical Information ---
if chemical_info:
    st.subheader("🧪 Chemical Information")
    st.write(chemical_info)

# --- Display Amino Acid Sequence ---
if sequence:
    st.subheader("🧬 Amino Acid Sequence (Extracted from PDB)")
    st.text_area("Sequence:", sequence, height=150)

# --- Display 3D Structure ---
if pdb_file and os.path.exists(pdb_file):
    st.subheader("🔬 Interactive 3D Structure Visualization")
    viewer = py3Dmol.view(width=1000, height=800)
    with open(pdb_file, "r") as f:
        pdb_content = f.read()
    viewer.addModel(pdb_content, "pdb")
    viewer.setStyle({"cartoon": {"color": "spectrum"}})
    viewer.zoomTo()
    st.components.v1.html(viewer._make_html(), height=800)
