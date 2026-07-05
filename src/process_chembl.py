import json
import os
from typing import Dict, List, Set, Tuple, Optional
import psycopg2
from psycopg2.extras import RealDictCursor


# =========================
# Why we extract ChEMBL drug-target interactions
# -----------------------------------
# Decagon requires drug-protein edges to build the (1,0) and (0,1) subgraphs
# in the multimodal network. ChEMBL provides curated drug-target data with:
#   - molecule_chembl_id (e.g., CHEMBL1234) as drug identifier
#   - UniProt accession (e.g., P04637) as protein/target identifier
#
# We filter for:
#   1. Human targets only (tax_id = 9606)
#   2. Direct binding mechanisms (high-confidence)
#   3. Single-protein targets (exclude protein complexes/families for now)
#
# Output format must match the UniProt ID space built from STRING PPI
# to enable consistent (0,0) protein-protein and (1,0) drug-protein edges.
# =========================

# Resolve project root (decagon-preprocessing/)
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_PROCESSED = os.path.join(PROJECT_ROOT, "data", "processed")


def connect_to_chembl(
    dbname: str = "chembl_36",
    user: str = "postgres",
    password: str = "",
    host: str = "localhost",
    port: int = 5432,
) -> psycopg2.extensions.connection:
    """
    Connect to local ChEMBL PostgreSQL database.

    Args:
        dbname: ChEMBL database name (default: chembl_36)
        user: PostgreSQL username
        password: PostgreSQL password
        host: PostgreSQL host
        port: PostgreSQL port

    Returns:
        psycopg2 connection object
    """
    conn = psycopg2.connect(
        dbname=dbname,
        user=user,
        password=password,
        host=host,
        port=port,
    )
    print(f"Connected to ChEMBL database: {dbname}@{host}:{port}")
    return conn


def extract_drug_target_interactions(
    conn: psycopg2.extensions.connection,
    action_types: Optional[List[str]] = None,
    max_rows: Optional[int] = None,
) -> List[Tuple[str, str, dict]]:
    """
    Extract drug-target interactions from ChEMBL.

    Query strategy:
      drug_mechanism  -> links drugs to targets + mechanism info
      molecule_dictionary -> get molecule_chembl_id
      target_dictionary -> get target info (filter single proteins)
      target_components -> link target to protein components
      component_sequences -> get UniProt accessions (filter human tax_id=9606)

    Args:
        conn: PostgreSQL connection
        action_types: Filter by mechanism action_types (e.g., ['INHIBITOR', 'ANTAGONIST'])
                      If None, keep all direct mechanisms.
        max_rows: Limit number of rows (for testing)

    Returns:
        List of (molecule_chembl_id, uniprot_accession, metadata_dict) tuples
    """

    # If action_types is provided, filter to those types.
    # If None (default), keep all mechanisms — maximizes target coverage
    # for intersection with STRING PPI.
    action_filter = ""
    if action_types:
        action_filter = "AND dm.action_type IN ({})".format(
            ",".join(f"'{at}'" for at in action_types)
        )

    limit_clause = f"LIMIT {max_rows}" if max_rows else ""

    query = f"""
    SELECT DISTINCT
        md.chembl_id AS molecule_chembl_id,
        cs.accession AS uniprot_accession,
        dm.action_type,
        dm.mechanism_of_action,
        td.target_type,
        cs.tax_id,
        td.pref_name AS target_name
    FROM
        drug_mechanism dm
    INNER JOIN
        molecule_dictionary md ON dm.molregno = md.molregno
    INNER JOIN
        target_dictionary td ON dm.tid = td.tid
    INNER JOIN
        target_components tc ON td.tid = tc.tid
    INNER JOIN
        component_sequences cs ON tc.component_id = cs.component_id
    WHERE
        cs.tax_id = 9606                          -- human only
        AND td.target_type = 'SINGLE PROTEIN'     -- exclude complexes/families
        AND cs.accession IS NOT NULL              -- must have UniProt accession
        AND dm.direct_interaction = 1              -- direct interactions only
        {action_filter}
    ORDER BY
        md.chembl_id, cs.accession
    {limit_clause}
    """

    print("Executing ChEMBL query...")
    print(f"  Action types filter: {action_types}")
    print(f"  Max rows: {max_rows if max_rows else 'unlimited'}")

    with conn.cursor(cursor_factory=RealDictCursor) as cur:
        cur.execute(query)
        rows = cur.fetchall()

    print(f"Fetched {len(rows):,} drug-target interaction rows")

    # Convert to list of tuples with metadata
    interactions = []
    for row in rows:
        molecule_id = row['molecule_chembl_id']
        uniprot = row['uniprot_accession']
        metadata = {
            'action_type': row['action_type'],
            'mechanism_of_action': row['mechanism_of_action'],
            'target_type': row['target_type'],
            'tax_id': row['tax_id'],
            'target_name': row['target_name'],
        }
        interactions.append((molecule_id, uniprot, metadata))

    return interactions


def extract_activity_based_targets(
    conn: psycopg2.extensions.connection,
    pchembl_min: float = 6.0,
    assay_types: Optional[List[str]] = None,
    standard_types: Optional[List[str]] = None,
    max_rows: Optional[int] = None,
) -> List[Tuple[str, str, dict]]:
    """
    Extract drug-target pairs from ChEMBL bioactivity data (activities + assays).

    This produces a broader "targets-all" set compared to the mechanism-based
    extraction. Filters for high-confidence binding/functional assay results.

    Args:
        conn: PostgreSQL connection
        pchembl_min: Minimum pChEMBL value threshold (default 6.0 = 1 uM)
        assay_types: Assay types to include (default: ['B', 'F'] for binding/functional)
        standard_types: Activity measurement types (default: Ki, IC50, EC50, Kd)
        max_rows: Limit number of rows (for testing)

    Returns:
        List of (molecule_chembl_id, uniprot_accession, metadata_dict) tuples
    """
    if assay_types is None:
        assay_types = ['B', 'F']
    if standard_types is None:
        standard_types = ['Ki', 'IC50', 'EC50', 'Kd']

    assay_filter = ",".join(f"'{at}'" for at in assay_types)
    std_filter = ",".join(f"'{st}'" for st in standard_types)
    limit_clause = f"LIMIT {max_rows}" if max_rows else ""

    # Aggregate per (drug, target): keep the best (max) pchembl_value
    query = f"""
    SELECT
        md.chembl_id AS molecule_chembl_id,
        cs.accession AS uniprot_accession,
        MAX(act.pchembl_value) AS best_pchembl,
        COUNT(*) AS n_measurements
    FROM
        activities act
    INNER JOIN
        molecule_dictionary md ON act.molregno = md.molregno
    INNER JOIN
        assays a ON act.assay_id = a.assay_id
    INNER JOIN
        target_dictionary td ON a.tid = td.tid
    INNER JOIN
        target_components tc ON td.tid = tc.tid
    INNER JOIN
        component_sequences cs ON tc.component_id = cs.component_id
    WHERE
        td.organism = 'Homo sapiens'
        AND td.target_type = 'SINGLE PROTEIN'
        AND cs.accession IS NOT NULL
        AND a.assay_type IN ({assay_filter})
        AND act.standard_type IN ({std_filter})
        AND act.standard_relation = '='
        AND act.pchembl_value >= {pchembl_min}
    GROUP BY
        md.chembl_id, cs.accession
    ORDER BY
        md.chembl_id, cs.accession
    {limit_clause}
    """

    print("Executing ChEMBL activities query...")
    print(f"  pChEMBL >= {pchembl_min}")
    print(f"  Assay types: {assay_types}")
    print(f"  Standard types: {standard_types}")
    print(f"  Max rows: {max_rows if max_rows else 'unlimited'}")

    with conn.cursor(cursor_factory=RealDictCursor) as cur:
        cur.execute(query)
        rows = cur.fetchall()

    print(f"Fetched {len(rows):,} activity-based drug-target pairs")

    interactions = []
    for row in rows:
        molecule_id = row['molecule_chembl_id']
        uniprot = row['uniprot_accession']
        metadata = {
            'source': 'activity',
            'best_pchembl': float(row['best_pchembl']),
            'n_measurements': int(row['n_measurements']),
        }
        interactions.append((molecule_id, uniprot, metadata))

    return interactions


def deduplicate_and_aggregate(
    interactions: List[Tuple[str, str, dict]]
) -> Tuple[List[Tuple[str, str]], Dict[Tuple[str, str], List[dict]]]:
    """
    Deduplicate drug-target pairs and aggregate metadata.

    Multiple rows may exist for the same (drug, target) pair if there are
    multiple mechanisms or evidence sources. We keep unique pairs and
    aggregate all metadata.

    Args:
        interactions: List of (molecule_id, uniprot, metadata) tuples

    Returns:
        unique_pairs: List of unique (molecule_id, uniprot) tuples
        pair_metadata: Dict mapping (molecule_id, uniprot) -> list of metadata dicts
    """
    pair_metadata: Dict[Tuple[str, str], List[dict]] = {}

    for molecule_id, uniprot, meta in interactions:
        key = (molecule_id, uniprot)
        if key not in pair_metadata:
            pair_metadata[key] = []
        pair_metadata[key].append(meta)

    unique_pairs = list(pair_metadata.keys())

    print(f"Unique (drug, target) pairs: {len(unique_pairs):,}")
    print(f"  Raw interaction rows: {len(interactions):,}")
    print(f"  Deduplication ratio: {len(interactions) / len(unique_pairs):.2f}x")

    return unique_pairs, pair_metadata


def analyze_coverage(
    unique_pairs: List[Tuple[str, str]],
    pair_metadata: Dict[Tuple[str, str], List[dict]],
) -> Dict:
    """
    Analyze drug and target coverage statistics.

    Returns summary statistics for thesis reporting.
    """
    drugs: Set[str] = set()
    targets: Set[str] = set()
    action_types_count: Dict[str, int] = {}

    for drug, target in unique_pairs:
        drugs.add(drug)
        targets.add(target)

        # Count action types / sources
        for meta in pair_metadata[(drug, target)]:
            action = meta.get('action_type') or meta.get('source', 'unknown')
            action_types_count[action] = action_types_count.get(action, 0) + 1

    stats = {
        'n_unique_drugs': len(drugs),
        'n_unique_targets': len(targets),
        'n_unique_pairs': len(unique_pairs),
        'action_type_distribution': action_types_count,
        'avg_targets_per_drug': len(unique_pairs) / len(drugs) if drugs else 0,
        'avg_drugs_per_target': len(unique_pairs) / len(targets) if targets else 0,
    }

    print("\n=== ChEMBL Drug-Target Coverage ===")
    print(f"  Unique drugs:   {stats['n_unique_drugs']:,}")
    print(f"  Unique targets: {stats['n_unique_targets']:,}")
    print(f"  Unique pairs:   {stats['n_unique_pairs']:,}")
    print(f"  Avg targets/drug: {stats['avg_targets_per_drug']:.2f}")
    print(f"  Avg drugs/target: {stats['avg_drugs_per_target']:.2f}")
    print("\n  Action type distribution:")
    for action, count in sorted(action_types_count.items(), key=lambda x: -x[1]):
        print(f"    {action:20s}: {count:6,} ({100*count/len(unique_pairs):.1f}%)")

    return stats


def build_id_mappings(
    unique_pairs: List[Tuple[str, str]]
) -> Tuple[Dict[str, int], Dict[str, int]]:
    """
    Build consistent integer index mappings for drugs and targets.

    Args:
        unique_pairs: List of (molecule_chembl_id, uniprot) tuples

    Returns:
        drug_to_idx: Dict mapping molecule_chembl_id -> integer index
        target_to_idx: Dict mapping uniprot_accession -> integer index
    """
    drug_to_idx: Dict[str, int] = {}
    target_to_idx: Dict[str, int] = {}

    for drug, target in unique_pairs:
        if drug not in drug_to_idx:
            drug_to_idx[drug] = len(drug_to_idx)
        if target not in target_to_idx:
            target_to_idx[target] = len(target_to_idx)

    print(f"\nBuilt ID mappings:")
    print(f"  Drugs:   {len(drug_to_idx):,} unique IDs")
    print(f"  Targets: {len(target_to_idx):,} unique IDs")

    return drug_to_idx, target_to_idx


def save_artifacts(
    out_dir: str,
    unique_pairs: List[Tuple[str, str]],
    drug_to_idx: Dict[str, int],
    target_to_idx: Dict[str, int],
    pair_metadata: Dict[Tuple[str, str], List[dict]],
    stats: Dict,
    config: Dict,
) -> None:
    """
    Save ChEMBL drug-target artifacts to disk.

    Saves:
      - drug_target_pairs.tsv: (molecule_chembl_id, uniprot_accession) pairs
      - drug_to_idx.json: Drug ID -> integer index mapping
      - target_to_idx.json: Target UniProt -> integer index mapping
      - pair_metadata.json: Full metadata for each pair
      - stats.json: Coverage statistics
      - meta.json: Processing configuration
    """
    os.makedirs(out_dir, exist_ok=True)

    # Save pairs as TSV (human-readable)
    pairs_path = os.path.join(out_dir, "drug_target_pairs.tsv")
    with open(pairs_path, "w", encoding="utf-8") as f:
        f.write("molecule_chembl_id\tuniprot_accession\n")
        for drug, target in sorted(unique_pairs):
            f.write(f"{drug}\t{target}\n")
    print(f"Saved: {pairs_path} ({len(unique_pairs):,} pairs)")

    # Save ID mappings
    drug_idx_path = os.path.join(out_dir, "drug_to_idx.json")
    with open(drug_idx_path, "w", encoding="utf-8") as f:
        json.dump(drug_to_idx, f)
    print(f"Saved: {drug_idx_path}")

    target_idx_path = os.path.join(out_dir, "target_to_idx.json")
    with open(target_idx_path, "w", encoding="utf-8") as f:
        json.dump(target_to_idx, f)
    print(f"Saved: {target_idx_path}")

    # Save metadata (convert tuple keys to strings for JSON)
    metadata_path = os.path.join(out_dir, "pair_metadata.json")
    serializable_metadata = {
        f"{drug}|{target}": metas
        for (drug, target), metas in pair_metadata.items()
    }
    with open(metadata_path, "w", encoding="utf-8") as f:
        json.dump(serializable_metadata, f, indent=2)
    print(f"Saved: {metadata_path}")

    # Save statistics
    stats_path = os.path.join(out_dir, "stats.json")
    with open(stats_path, "w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2)
    print(f"Saved: {stats_path}")

    # Save configuration
    meta_path = os.path.join(out_dir, "meta.json")
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(config, f, indent=2)
    print(f"Saved: {meta_path}")


if __name__ == "__main__":
    # =====================================
    # Configuration
    # =====================================
    import sys as _sys
    _src = os.path.dirname(os.path.abspath(__file__))
    if _src not in _sys.path:
        _sys.path.insert(0, _src)
    import config as _cfg

    # PostgreSQL connection — all values configurable via environment variables
    DB_CONFIG = {
        "dbname":   _cfg.CHEMBL_DB_NAME,
        "user":     _cfg.CHEMBL_DB_USER,
        "password": _cfg.CHEMBL_DB_PASSWORD,
        "host":     _cfg.CHEMBL_DB_HOST,
        "port":     _cfg.CHEMBL_DB_PORT,
    }

    # Action types: None = keep all mechanisms (maximizes target coverage)
    ACTION_TYPES = None

    # Activities-based extraction config
    PCHEMBL_MIN = 6.0  # 1 uM threshold

    # For initial testing, limit rows; set to None for full extraction
    MAX_ROWS = None

    # Output directories
    OUT_DIR_STRICT = os.path.join(DATA_PROCESSED, "chembl")
    OUT_DIR_ALL = os.path.join(DATA_PROCESSED, "chembl_all")

    # =====================================
    # Execution
    # =====================================

    print("=" * 60)
    print("ChEMBL Drug-Target Extraction")
    print("=" * 60)

    # 1. Connect to database
    conn = connect_to_chembl(**DB_CONFIG)

    try:
        # ── Strict targets (mechanism-based) ──────────────────────
        print("\n--- Strict targets (drug_mechanism) ---")
        mechanism_interactions = extract_drug_target_interactions(
            conn,
            action_types=ACTION_TYPES,
            max_rows=MAX_ROWS,
        )

        mech_pairs, mech_metadata = deduplicate_and_aggregate(mechanism_interactions)
        mech_stats = analyze_coverage(mech_pairs, mech_metadata)
        mech_drug_idx, mech_target_idx = build_id_mappings(mech_pairs)

        config_strict = {
            "db_config": {k: v for k, v in DB_CONFIG.items() if k != "password"},
            "source": "drug_mechanism",
            "action_types": ACTION_TYPES,
            "direct_interaction": True,
            "filters": {"tax_id": 9606, "target_type": "SINGLE PROTEIN"},
        }

        save_artifacts(
            OUT_DIR_STRICT, mech_pairs, mech_drug_idx, mech_target_idx,
            mech_metadata, mech_stats, config_strict,
        )

        # ── Broad targets (activities-based) ──────────────────────
        print("\n--- Broad targets (activities) ---")
        activity_interactions = extract_activity_based_targets(
            conn,
            pchembl_min=PCHEMBL_MIN,
            max_rows=MAX_ROWS,
        )

        # ── Union: mechanism + activities ─────────────────────────
        print("\n--- Union (mechanism + activities) ---")
        all_interactions = mechanism_interactions + activity_interactions
        all_pairs, all_metadata = deduplicate_and_aggregate(all_interactions)
        all_stats = analyze_coverage(all_pairs, all_metadata)
        all_drug_idx, all_target_idx = build_id_mappings(all_pairs)

        # Compare
        mech_pair_set = set(mech_pairs)
        all_pair_set = set(all_pairs)
        activity_only = all_pair_set - mech_pair_set
        mech_targets = {t for _, t in mech_pairs}
        all_targets = {t for _, t in all_pairs}

        print(f"\n=== Coverage comparison ===")
        print(f"  Mechanism-only pairs:  {len(mech_pair_set):,}")
        print(f"  Union pairs:           {len(all_pair_set):,}")
        print(f"  Added by activities:   {len(activity_only):,}")
        print(f"  Mechanism targets:     {len(mech_targets):,}")
        print(f"  Union targets:         {len(all_targets):,}")
        print(f"  New targets from activities: {len(all_targets - mech_targets):,}")

        config_all = {
            "db_config": {k: v for k, v in DB_CONFIG.items() if k != "password"},
            "source": "union(drug_mechanism, activities)",
            "mechanism_action_types": ACTION_TYPES,
            "mechanism_direct_interaction": True,
            "activity_pchembl_min": PCHEMBL_MIN,
            "activity_assay_types": ["B", "F"],
            "activity_standard_types": ["Ki", "IC50", "EC50", "Kd"],
            "filters": {"tax_id": 9606, "target_type": "SINGLE PROTEIN"},
        }

        save_artifacts(
            OUT_DIR_ALL, all_pairs, all_drug_idx, all_target_idx,
            all_metadata, all_stats, config_all,
        )

        print("\n" + "=" * 60)
        print("ChEMBL processing complete!")
        print("=" * 60)

    finally:
        conn.close()
        print("Database connection closed.")
