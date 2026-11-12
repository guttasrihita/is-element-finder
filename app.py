import io
import pandas as pd
import streamlit as st
from Bio import SeqIO
from algorithm import clean_sequence, find_candidates, analyse_sequence

st.set_page_config(page_title="IS Element Finder", page_icon="🧬", layout="wide")
st.title("🧬 IS Element Finder")
st.caption("Screen DNA sequences for candidate bacterial insertion sequence (IS) elements.")

with st.sidebar:
    st.header("Analysis settings")
    min_len = st.number_input("Minimum candidate length (bp)", 40, 10000, 700, 10)
    max_len = st.number_input("Maximum candidate length (bp)", int(min_len), 50000, max(2500, int(min_len)), 10)
    tir_len = st.number_input("TIR length (bp)", 1, 100, 20, 1)
    mismatches = st.number_input("Maximum TIR mismatches", 0, int(tir_len), min(2, int(tir_len)), 1)
    tsd_len = st.number_input("TSD length to check (bp)", 1, 30, 5, 1)

st.markdown("Paste a DNA sequence or upload a FASTA file. The tool checks terminal inverted repeat (TIR) similarity and reports whether a matching target-site duplication (TSD) is detected.")
tab1, tab2 = st.tabs(["Paste DNA sequence", "Upload FASTA"])

records = []
with tab1:
    sequence_text = st.text_area("DNA sequence (A, C, G, T, optionally N)", height=180, placeholder="Paste DNA sequence here...")
    record_name = st.text_input("Sequence name", "Pasted_sequence")
    if st.button("Analyse pasted sequence", type="primary"):
        if not sequence_text.strip():
            st.warning("Paste a DNA sequence first.")
        else:
            try:
                records = [(record_name.strip() or "Pasted_sequence", clean_sequence(sequence_text))]
            except ValueError as e:
                st.error(str(e))

with tab2:
    uploaded = st.file_uploader("Upload FASTA file", type=["fasta", "fa", "fna", "txt"])
    if st.button("Analyse FASTA file", type="primary"):
        if uploaded is None:
            st.warning("Upload a FASTA file first.")
        else:
            try:
                parsed = list(SeqIO.parse(io.StringIO(uploaded.getvalue().decode("utf-8")), "fasta"))
                if not parsed:
                    st.error("No FASTA records found. Check the file format.")
                else:
                    records = [(r.id, clean_sequence(str(r.seq))) for r in parsed]
            except (UnicodeDecodeError, ValueError) as e:
                st.error(f"Could not read FASTA file: {e}")

if records:
    rows = []
    st.subheader("Analysis results")
    for record_id, sequence in records:
        try:
            candidates = find_candidates(sequence, int(min_len), int(max_len), int(tir_len), int(mismatches))
            results = analyse_sequence(sequence, candidates, int(tsd_len))
            st.markdown(f"**{record_id}** · {len(sequence):,} bp")
            if not results:
                st.info("No candidate IS elements detected with these settings.")
            else:
                st.success(f"{len(results)} candidate(s) detected.")
                for number, hit in enumerate(results, 1):
                    with st.expander(f"Candidate {number}: {hit['start']}–{hit['end']} ({hit['length']} bp)", expanded=True):
                        c1, c2, c3 = st.columns(3)
                        c1.metric("Element length", f"{hit['length']:,} bp")
                        c2.metric("TIR mismatches", hit["tir_mismatches"])
                        c3.metric("TSD", hit["target_site_duplication"])
                        st.write("Left TIR:", hit["left_tir"])
                        st.write("Right TIR:", hit["right_tir"])
                rows.extend({"record_id": record_id, **hit} for hit in results)
        except ValueError as e:
            st.error(f"{record_id}: {e}")
    if rows:
        df = pd.DataFrame(rows)
        st.dataframe(df, use_container_width=True, hide_index=True)
        st.download_button("Download results as CSV", df.to_csv(index=False).encode("utf-8"),
                           "is_element_finder_results.csv", "text/csv")

st.divider()
st.caption("Educational research prototype. These are heuristic candidate detections, not confirmed IS annotations. Synthetic examples test software behaviour but do not establish accuracy on real bacterial genomes.")
