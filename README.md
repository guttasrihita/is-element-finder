# IS Element Finder

## About the Project

IS Element Finder is a Python-based web application designed to identify potential bacterial Insertion Sequence (IS) elements from DNA sequences.

The application analyses DNA sequences by looking for Terminal Inverted Repeats (TIRs) and checking for possible Target Site Duplications (TSDs).

## Features

- Enter a DNA sequence directly into the application.
- Upload DNA sequences in FASTA format.
- Identify potential IS elements based on TIR similarity.
- Configure sequence length and mismatch parameters.
- View candidate positions and sequence details.
- Check for possible target site duplications.
- Download analysis results as a CSV file.

## Technologies Used

- Python
- Biopython
- Streamlit
- Pandas

## How to Use

1. Open the IS Element Finder web application.
2. Paste a DNA sequence or upload a FASTA file.
3. Select the required analysis settings.
4. Click the Analyse button.
5. View the detected candidate regions and their details.
6. Download the results if required.

## Project Structure

- `app.py` – Web application interface.
- `algorithm.py` – DNA sequence analysis and IS candidate detection.
- `requirements.txt` – Required Python libraries.
- `synthetic_demo.fasta` – Sample sequence for testing.

## Note

This application identifies potential IS elements using a heuristic approach based on terminal inverted repeats and target site duplications. The detected candidates should not be considered experimentally confirmed IS elements. Further validation using real bacterial sequences and reference databases is required.

## Project Team

- Shraddha Padmaraju
- Srihita Gutta
