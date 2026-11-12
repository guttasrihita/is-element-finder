import argparse
import csv
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq


def clean_sequence(sequence):
    sequence = "".join(sequence.split()).upper()
    invalid = sorted(set(sequence) - set("ACGTN"))
    if invalid:
        raise ValueError("Invalid DNA characters: " + ", ".join(invalid))
    if not sequence:
        raise ValueError("The DNA sequence is empty.")
    return sequence


def read_fasta(filename):
    path = Path(filename)
    if not path.is_file():
        raise FileNotFoundError(f"FASTA file not found: {filename}")

    records = []
    for record in SeqIO.parse(str(path), "fasta"):
        records.append((record.id, clean_sequence(str(record.seq))))

    if not records:
        raise ValueError("No FASTA records found.")
    return records


def reverse_complement(sequence):
    return str(Seq(sequence).reverse_complement())


def hamming_distance(first, second):
    if len(first) != len(second):
        raise ValueError("Sequences must have equal lengths.")
    return sum(a != b for a, b in zip(first, second))


def find_candidates(sequence, min_len=700, max_len=2500,
                    tir_len=20, max_mismatches=2):
    sequence = clean_sequence(sequence)

    if min_len <= 0 or max_len < min_len:
        raise ValueError("Use 0 < minimum length <= maximum length.")
    if tir_len <= 0 or 2 * tir_len > min_len:
        raise ValueError("TIR length must be positive and fit within min_len.")
    if max_mismatches < 0 or max_mismatches > tir_len:
        raise ValueError("Mismatch limit must be between 0 and TIR length.")

    hits = []
    n = len(sequence)

    for start0 in range(max(0, n - min_len + 1)):
        left_tir = sequence[start0:start0 + tir_len]
        if len(left_tir) != tir_len or "N" in left_tir:
            continue

        longest = min(max_len, n - start0)
        for element_length in range(min_len, longest + 1):
            end0 = start0 + element_length
            right_tir = sequence[end0 - tir_len:end0]

            if len(right_tir) != tir_len or "N" in right_tir:
                continue

            mismatches = hamming_distance(
                reverse_complement(left_tir), right_tir
            )

            if mismatches <= max_mismatches:
                hits.append({
                    "start": start0 + 1,
                    "end": end0,
                    "length": element_length,
                    "left_tir": left_tir,
                    "right_tir": right_tir,
                    "tir_mismatches": mismatches
                })
                break

    hits.sort(key=lambda hit: (
        hit["tir_mismatches"], -hit["length"], hit["start"]
    ))

    selected = []
    for hit in hits:
        overlaps = any(
            not (hit["end"] < chosen["start"] or
                 hit["start"] > chosen["end"])
            for chosen in selected
        )
        if not overlaps:
            selected.append(hit)

    return sorted(selected, key=lambda hit: hit["start"])


def check_target_site_duplication(sequence, start, end, duplication_len=5):
    sequence = clean_sequence(sequence)

    if duplication_len <= 0:
        raise ValueError("Duplication length must be positive.")

    left_start = start - 1 - duplication_len
    left_end = start - 1
    right_start = end
    right_end = end + duplication_len

    if left_start < 0 or right_end > len(sequence):
        return None

    left_repeat = sequence[left_start:left_end]
    right_repeat = sequence[right_start:right_end]

    if "N" not in left_repeat and left_repeat == right_repeat:
        return left_repeat
    return None


def analyse_sequence(sequence, candidates, duplication_len=5):
    results = []
    for hit in candidates:
        row = dict(hit)
        repeat = check_target_site_duplication(
            sequence, hit["start"], hit["end"], duplication_len
        )
        row["target_site_duplication"] = repeat or "Not detected"
        results.append(row)
    return results


def print_results(record_id, sequence, results):
    print(f"\\nRecord: {record_id}")
    print(f"Sequence length: {len(sequence)} bp")

    if not results:
        print("No candidate IS elements detected.")
        return

    print(f"Candidate IS elements: {len(results)}")
    print("-" * 86)
    print(f"{'No.':<5}{'Start':<10}{'End':<10}{'Length':<10}"
          f"{'TIR mismatches':<17}{'TSD':<15}")
    print("-" * 86)

    for index, hit in enumerate(results, 1):
        print(f"{index:<5}{hit['start']:<10}{hit['end']:<10}"
              f"{hit['length']:<10}{hit['tir_mismatches']:<17}"
              f"{hit['target_site_duplication']:<15}")
        print("  Left TIR :", hit["left_tir"])
        print("  Right TIR:", hit["right_tir"])


def save_csv(all_results, filename):
    fields = [
        "record_id", "start", "end", "length", "left_tir",
        "right_tir", "tir_mismatches", "target_site_duplication"
    ]
    with open(filename, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(all_results)


def main():
    parser = argparse.ArgumentParser(
        description="Screen DNA sequences for candidate bacterial IS elements."
    )
    inputs = parser.add_mutually_exclusive_group()
    inputs.add_argument("-f", "--fasta", help="Input FASTA file.")
    inputs.add_argument("-s", "--sequence", help="DNA sequence as text.")
    parser.add_argument("--min-len", type=int, default=700)
    parser.add_argument("--max-len", type=int, default=2500)
    parser.add_argument("--tir-len", type=int, default=20)
    parser.add_argument("--mismatches", type=int, default=2)
    parser.add_argument("--tsd-len", type=int, default=5)
    parser.add_argument("--csv", help="Optional CSV output file.")
    args = parser.parse_args()

    try:
        if args.fasta:
            records = read_fasta(args.fasta)
        elif args.sequence:
            records = [("input_sequence", clean_sequence(args.sequence))]
        else:
            records = [("input_sequence",
                        clean_sequence(input("Enter DNA sequence: ")))]

        all_results = []
        for record_id, sequence in records:
            candidates = find_candidates(
                sequence, args.min_len, args.max_len,
                args.tir_len, args.mismatches
            )
            results = analyse_sequence(sequence, candidates, args.tsd_len)
            print_results(record_id, sequence, results)
            all_results.extend(
                {"record_id": record_id, **result} for result in results
            )

        if args.csv:
            save_csv(all_results, args.csv)
            print(f"\\nResults saved to: {args.csv}")

    except (ValueError, FileNotFoundError) as error:
        raise SystemExit(f"Error: {error}")


if __name__ == "__main__":
    main()
