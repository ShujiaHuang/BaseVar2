"""
Configuration for synthetic test data generation.

Defines chromosomes, ground-truth variants, sample groups,
and simulation parameters for basevar caller test data.
"""

# ---------------------------------------------------------------------------
# Random seed (fixed for reproducibility)
# ---------------------------------------------------------------------------
DEFAULT_SEED = 42

# ---------------------------------------------------------------------------
# Reference genome design
# ---------------------------------------------------------------------------
CHROMOSOMES = {
    "chr1": 2000,  # Main test chromosome, carries most variants
    "chr2": 1500,  # Multi-chromosome region test
    "chrX": 1000,  # Sex chromosome test
    "chrY":  500,  # Reserved (not used for caller tests)
}
GC_CONTENT = 0.45  # Target GC content ~45%

# ---------------------------------------------------------------------------
# Ground truth variants (11 variants)
#
# IMPORTANT: The REF/ALT values below are DESIGN-TIME placeholders.
# During data generation, `adapt_variants_to_ref()` in reference.py
# automatically adjusts REF/ALT to match the actual random reference
# sequence. The adapted values are written to ground_truth_variants.tsv.
#
# For multi-allelic sites (v5), ALT is comma-separated.
# AF values are per-group allele frequencies.
# For multi-allelic v5: AF = (AF_alt1, AF_alt2) per group.
# ---------------------------------------------------------------------------
VARIANTS = [
    # ID    CHROM   POS(1-based)  REF   ALT          TYPE            GroupA_AF  GroupB_AF  GroupC_AF
    {"id": "v1",  "chrom": "chr1", "pos": 200,  "ref": "A",   "alt": "G",       "type": "SNP",   "af": {"GroupA": 0.5, "GroupB": 0.5, "GroupC": 0.5}},
    {"id": "v2",  "chrom": "chr1", "pos": 400,  "ref": "C",   "alt": "T",       "type": "SNP",   "af": {"GroupA": 0.1, "GroupB": 0.5, "GroupC": 0.8}},
    {"id": "v3",  "chrom": "chr1", "pos": 600,  "ref": "G",   "alt": "A",       "type": "SNP",   "af": {"GroupA": 0.0, "GroupB": 0.3, "GroupC": 0.0}},
    {"id": "v4",  "chrom": "chr1", "pos": 800,  "ref": "T",   "alt": "C",       "type": "SNP",   "af": {"GroupA": 0.3, "GroupB": 0.0, "GroupC": 0.6}},
    {"id": "v5",  "chrom": "chr1", "pos": 1000, "ref": "A",   "alt": "G,T",     "type": "multi-allelic", "af": {"GroupA": (0.3, 0.1), "GroupB": (0.2, 0.05), "GroupC": (0.4, 0.1)}},
    {"id": "v6",  "chrom": "chr1", "pos": 1200, "ref": "GA",  "alt": "G",       "type": "del",   "af": {"GroupA": 0.2, "GroupB": 0.4, "GroupC": 0.1}},
    {"id": "v7",  "chrom": "chr1", "pos": 1400, "ref": "C",   "alt": "CA",      "type": "ins",   "af": {"GroupA": 0.2, "GroupB": 0.2, "GroupC": 0.5}},
    {"id": "v8",  "chrom": "chr1", "pos": 1600, "ref": "GAT", "alt": "G",       "type": "del",   "af": {"GroupA": 0.1, "GroupB": 0.3, "GroupC": 0.2}},
    {"id": "v9",  "chrom": "chr2", "pos": 300,  "ref": "A",   "alt": "C",       "type": "SNP",   "af": {"GroupA": 0.4, "GroupB": 0.4, "GroupC": 0.4}},
    {"id": "v10", "chrom": "chr2", "pos": 700,  "ref": "T",   "alt": "G",       "type": "SNP",   "af": {"GroupA": 0.2, "GroupB": 0.6, "GroupC": 0.0}},
    {"id": "v11", "chrom": "chrX", "pos": 200,  "ref": "G",   "alt": "A",       "type": "SNP",   "af": {"GroupA": 0.3, "GroupB": 0.5, "GroupC": 0.2}},
]

# ---------------------------------------------------------------------------
# Sample design (50 samples, 3 groups)
# ---------------------------------------------------------------------------
SAMPLE_GROUPS = {
    "GroupA": [f"sampleA{i:02d}" for i in range(1, 21)],   # 20 samples
    "GroupB": [f"sampleB{i:02d}" for i in range(1, 21)],   # 20 samples
    "GroupC": [f"sampleC{i:02d}" for i in range(1, 11)],   # 10 samples
}

# Samples with reduced coverage (~1x instead of ~8x)
LOW_COVERAGE_SAMPLES = ["sampleC01", "sampleC02", "sampleC03", "sampleC04"]

# Samples with high coverage (~30x) for edge case testing
HIGH_COVERAGE_SAMPLES = ["sampleA01"]
HIGH_COVERAGE = 30

# ---------------------------------------------------------------------------
# Read simulation parameters
# ---------------------------------------------------------------------------
READ_LENGTH = 100        # Single-end reads
NORMAL_COVERAGE = 8      # Target coverage for most samples
LOW_COVERAGE = 1         # Target coverage for LOW_COVERAGE_SAMPLES
MAPQ_RANGE = (20, 60)    # Uniform distribution for mapping quality
BASEQ_MEAN = 30          # Mean base quality
BASEQ_SD = 5             # Standard deviation for base quality
BASEQ_CLIP = (10, 40)    # Clip base quality to this range
FRAGMENT_WINDOW = 200    # Variant reads generated within +/- this many bp of variant

# ---------------------------------------------------------------------------
# Derived helpers
# ---------------------------------------------------------------------------
def get_all_samples():
    """Return flat list of all sample IDs."""
    samples = []
    for group_samples in SAMPLE_GROUPS.values():
        samples.extend(group_samples)
    return samples

def get_sample_group(sample_id):
    """Return the group name for a given sample ID."""
    for group, samples in SAMPLE_GROUPS.items():
        if sample_id in samples:
            return group
    raise ValueError(f"Unknown sample: {sample_id}")

def get_target_coverage(sample_id):
    """Return target coverage for a sample."""
    if sample_id in LOW_COVERAGE_SAMPLES:
        return LOW_COVERAGE
    if sample_id in HIGH_COVERAGE_SAMPLES:
        return HIGH_COVERAGE
    return NORMAL_COVERAGE

def get_variant_af(variant, group):
    """
    Return allele frequency for a variant in a given group.
    For multi-allelic variants, returns a tuple of AFs.
    """
    af = variant["af"][group]
    return af
