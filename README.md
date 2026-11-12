# IS Element Finder

A Streamlit web application for screening DNA sequences for candidate bacterial insertion sequence (IS) elements.

## Features
- Paste a DNA sequence or upload a FASTA file.
- Configure candidate length, TIR length, mismatch limit and TSD length.
- View candidate coordinates, length, TIR sequences, mismatch count and TSD result.
- Download results as CSV.
- Supports multiple FASTA records.

## Scientific limitation
This is a heuristic candidate finder, not a validated IS annotation pipeline. TIR similarity and a possible TSD do not confirm an IS element. Synthetic test sequences are for controlled software testing only; validation on real bacterial datasets and curated IS references is future work.

## Run locally
Install Python 3.10+, then run:

```bash
pip install -r requirements.txt
streamlit run app.py
```

Open the local URL shown in the terminal (usually http://localhost:8501).

## Publish to GitHub
1. Create a **public** GitHub repository named `is-element-finder`.
2. Upload `app.py`, `algorithm.py`, `requirements.txt`, and `README.md` to the repository root.
3. Commit the files.

## Deploy a public live URL
1. Sign in to Streamlit Community Cloud using your GitHub account.
2. Select **Create app** / **Deploy an app**.
3. Choose repository `is-element-finder`, branch `main`, and main file `app.py`.
4. Deploy and wait for the build to complete.
5. Copy the generated `https://....streamlit.app` URL and test it in an incognito/private browser window.

The GitHub repository URL and live app URL are different. The live URL is only created after successful deployment.
