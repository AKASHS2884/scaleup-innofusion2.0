import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw, RDKFingerprint

from utils.dft import predict_dft_properties
from utils.admet import check_admet_properties
from utils.molecule_ops import generate_permutations, generate_combinations, apply_mutation
from utils.genetics import simulate_compatibility, genes, structures

#  App Configuration
st.set_page_config(page_title="Quantum Drug Innovator - ScaleUp", layout="wide")
st.title("🔬 QUANTUM DRUG INNOVATOR - SCALE UP")

#  SMILES Input
smiles = st.text_input("Enter SMILES HERE")

# 🧪 Main Logic
if smiles:
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        st.image(Draw.MolToImage(mol, size=(300, 300)))
        st.success("✅ Molecule successfully read!")

        # 🧬 Molecular Fingerprint
        st.subheader("🧬 RDKit Fingerprint")
        fingerprint = RDKFingerprint(mol)
        on_bits = list(fingerprint.GetOnBits())
        st.code(f"ON bits: {on_bits}")

        # 🧪 ADMET Evaluation
        st.subheader("🧪 ADMET Evaluation")
        admet_result = check_admet_properties(smiles)
        st.json(admet_result)

        if admet_result.get("Passes ADMET"):
            st.success("✅ Molecule passes ADMET properties!")
        else:
            st.warning("⚠️ Molecule may not meet ADMET criteria.")

        # 🧬 Molecule Variants & Evaluation
        st.subheader("🧬 Molecule Variants & ADMET Comparison")
        with st.spinner("Generating and evaluating variants..."):
            variant_list = []

            # Generate variants
            variant_list.extend(generate_permutations(smiles))
            variant_list.extend(generate_combinations(smiles, ["C", "CC", "CO"]))
            variant_list.append(apply_mutation(smiles))

            for idx, variant in enumerate(variant_list):
                st.markdown(f"---\n### 🔹 Variant {idx + 1}")
                st.write(f"**SMILES:** `{variant}`")

                variant_mol = Chem.MolFromSmiles(variant)
                if variant_mol:
                    st.image(Draw.MolToImage(variant_mol, size=(250, 250)))
                    variant_fp = RDKFingerprint(variant_mol)
                    st.code(f"ON bits: {list(variant_fp.GetOnBits())}")

                admet_result = check_admet_properties(variant)
                st.json(admet_result)

        # 🧬 Gene-DNA Merge Compatibility (Bonus Simulator)
        st.subheader("🧬 Gene-DNA Merge Compatibility Simulator")

        lang = st.selectbox("🌍 Verdict Language", ["English", "தமிழ்"])
        selected_gene = st.selectbox("🔬 Choose a Gene", genes)
        selected_structure = st.selectbox("🧪 Choose a Human RNA/DNA Structure", structures)

        if st.button("Analyze Merge Potential"):
            score, verdict = simulate_compatibility(selected_gene, selected_structure)
            st.metric(label="Compatibility Score", value=score)

            verdict_text = {
                "High Merge Potential": "Merge பண்ணலாம்! 🔥",
                "Low Merge Potential": "அவ்வளவாக ஓட்டமில்லை 😅"
            }.get(verdict, verdict) if lang == "தமிழ்" else verdict

            st.success(f"Verdict: {verdict_text}")
    else:
        st.error("❌ Invalid SMILES notation. Please try again.")
else:
    st.info("ℹ️ Waiting for SMILES input to connect to OpenFrame.")