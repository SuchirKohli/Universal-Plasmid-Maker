from collections import Counter
import numpy as np
from scipy.signal import find_peaks


GC_WINDOW = 500
KMER = 9
ENRICH_WINDOW = 500
DIST_TOL = 1000
ORI_WINDOW = 300   # final ORI extraction size


def gc_skew(seq, window):
    sk = []
    for i in range(0, len(seq) - window + 1, window):
        chunk = seq[i:i+window]
        g, c = chunk.count("G"), chunk.count("C")
        sk.append((g - c) / (g + c) if g + c else 0)
    return np.array(sk)


def kmer_counts(seq, k):
    return Counter(seq[i:i+k] for i in range(len(seq) - k + 1))


def enrichment_profile(seq, bg_freq, k, window):
    scores = []
    for i in range(0, len(seq) - window + 1, window):
        chunk = seq[i:i+window]
        local_count = kmer_counts(chunk, k)
        w_total = sum(local_count.values())
        score = 0
        for mer, cnt in local_count.items():
            if mer in bg_freq and bg_freq[mer] > 0:
                score += (cnt / w_total) / bg_freq[mer]
        scores.append(score)
    return np.array(scores)


def predict_ori(sequence: str) -> dict:
    """
    Predict ORI using GC skew + k-mer enrichment.
    Returns a bounded ORI region.
    """

    seq = sequence.upper()
    seq_len = len(seq)

    # GC skew
    skew = gc_skew(seq, GC_WINDOW)

    # k-mer background
    bg_counts = kmer_counts(seq, KMER)
    bg_total = sum(bg_counts.values())
    bg_freq = {k: v / bg_total for k, v in bg_counts.items()}

    # enrichment
    enrichment = enrichment_profile(seq, bg_freq, KMER, ENRICH_WINDOW)

    # peaks in enrichment
    peaks, _ = find_peaks(enrichment, prominence=1)

    # GC skew sign change
    sign_change = np.where(np.diff(np.sign(skew)) != 0)[0]

    # candidate ORI positions
    candidates = []
    for p in peaks:
        coord = p * ENRICH_WINDOW
        if any(abs(coord - sc * GC_WINDOW) < DIST_TOL for sc in sign_change):
            candidates.append(coord)

    # fallback if no consensus
    if not candidates:
        # AT-rich fallback (important for small plasmids)
        best_start = 0
        best_score = -1
        for i in range(0, seq_len - ORI_WINDOW):
            frag = seq[i:i+ORI_WINDOW]
            score = frag.count("A") + frag.count("T")
            if score > best_score:
                best_score = score
                best_start = i
        ori_start = best_start
    else:
        ori_start = candidates[0]

    ori_end = min(ori_start + ORI_WINDOW, seq_len)

    return {
        "ori_start": ori_start,
        "ori_end": ori_end,
        "ori_sequence": seq[ori_start:ori_end],
        "method": "GC skew + k-mer enrichment"
    }

def predict_ori_candidates(sequence: str, top_n: int = 5):
    """
    Return top N ORI candidate regions with enriched k-mers.
    Intended for exploratory analysis and notebooks.

    Returns
    -------
    list of dicts, each with:
        - start
        - end
        - sequence
        - enriched_kmers (list of (kmer, score))
    """

    seq = sequence.upper()
    seq_len = len(seq)

    # --- reuse internal logic ---
    skew = gc_skew(seq, GC_WINDOW)

    bg_counts = kmer_counts(seq, KMER)
    bg_total = sum(bg_counts.values())
    bg_freq = {k: v / bg_total for k, v in bg_counts.items()}

    enrichment = enrichment_profile(seq, bg_freq, KMER, ENRICH_WINDOW)
    peaks, properties = find_peaks(enrichment, prominence=1)

    # GC skew sign changes
    sign_change = np.where(np.diff(np.sign(skew)) != 0)[0]

    candidates = []

    for p in peaks:
        coord = p * ENRICH_WINDOW

        if any(abs(coord - sc * GC_WINDOW) < DIST_TOL for sc in sign_change):
            start = coord
            end = min(start + ORI_WINDOW, seq_len)
            window_seq = seq[start:end]

            # compute enriched kmers in this window
            local_counts = kmer_counts(window_seq, KMER)
            local_total = sum(local_counts.values())

            kmer_scores = []
            for mer, cnt in local_counts.items():
                if mer in bg_freq and bg_freq[mer] > 0:
                    score = (cnt / local_total) / bg_freq[mer]
                    kmer_scores.append((mer, score))

            kmer_scores.sort(key=lambda x: x[1], reverse=True)

            candidates.append({
                "start": start,
                "end": end,
                "sequence": window_seq,
                "enriched_kmers": kmer_scores[:5]
            })

    # Fallback: AT-rich windows if no candidates
    if not candidates:
        for i in range(0, seq_len - ORI_WINDOW, ORI_WINDOW):
            window_seq = seq[i:i+ORI_WINDOW]
            at_score = window_seq.count("A") + window_seq.count("T")

            candidates.append({
                "start": i,
                "end": i + ORI_WINDOW,
                "sequence": window_seq,
                "enriched_kmers": [],
                "at_score": at_score
            })

        candidates.sort(key=lambda x: x["at_score"], reverse=True)

    return candidates[:top_n]
