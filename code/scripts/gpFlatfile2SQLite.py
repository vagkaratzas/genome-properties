#!/usr/bin/env python3
"""
Parse a Genome Properties flatfile (genomeProperties.txt) and export
its contents into a SQLite database.

The schema mirrors the original MySQL database used by gpFlatfile2DB.pl,
with 10 tables covering genome properties, steps, evidence (InterPro and
GenProp), GO term mappings, literature references, and database links.

GO term labels/categories are fetched from the EBI QuickGO REST API on
first encounter and cached in the database for subsequent runs.

Usage:
    python gpFlatfile2SQLite.py flatfiles/genomeProperties.txt -o genome_properties.db
"""

import argparse
import os
import re
import sqlite3
import sys
import time

import requests

# ---------------------------------------------------------------------------
# Schema  (mirrors the original MySQL / DBIx::Class schema)
# ---------------------------------------------------------------------------

SCHEMA_SQL = """\
CREATE TABLE IF NOT EXISTS genome_property (
    accession   TEXT    NOT NULL PRIMARY KEY,
    description TEXT    NOT NULL,
    type        TEXT    NOT NULL,
    author      TEXT,
    threshold   INTEGER NOT NULL DEFAULT 0,
    comment     TEXT,
    private     TEXT,
    ispublic    INTEGER NOT NULL DEFAULT 0,
    checked     INTEGER NOT NULL DEFAULT 0
);

CREATE TABLE IF NOT EXISTS gp_step (
    auto_step          INTEGER PRIMARY KEY AUTOINCREMENT,
    gp_accession       TEXT    NOT NULL REFERENCES genome_property(accession) ON DELETE CASCADE,
    step_number        INTEGER NOT NULL,
    step_id            TEXT    NOT NULL,
    step_display_name  TEXT,
    required           INTEGER NOT NULL
);

CREATE TABLE IF NOT EXISTS gp_step_evidence_ipr (
    auto_ipr_step  INTEGER PRIMARY KEY AUTOINCREMENT,
    auto_step      INTEGER NOT NULL REFERENCES gp_step(auto_step) ON DELETE CASCADE,
    interpro_acc   TEXT    NOT NULL,
    signature_acc  TEXT    NOT NULL,
    sufficient     INTEGER NOT NULL
);

CREATE TABLE IF NOT EXISTS gp_step_evidence_gp (
    auto_gp_step   INTEGER PRIMARY KEY AUTOINCREMENT,
    auto_step      INTEGER NOT NULL REFERENCES gp_step(auto_step) ON DELETE CASCADE,
    gp_accession   TEXT    NOT NULL REFERENCES genome_property(accession) ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS go_terms (
    go_id    TEXT NOT NULL PRIMARY KEY,
    term     TEXT NOT NULL,
    category TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS gp_step_to_go (
    auto_gp_step INTEGER NOT NULL REFERENCES gp_step_evidence_gp(auto_gp_step) ON DELETE CASCADE,
    go_id        TEXT    NOT NULL REFERENCES go_terms(go_id) ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS ipr_step_to_go (
    auto_ipr_step INTEGER NOT NULL REFERENCES gp_step_evidence_ipr(auto_ipr_step) ON DELETE CASCADE,
    go_id         TEXT    NOT NULL REFERENCES go_terms(go_id) ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS literature_reference (
    pmid    INTEGER NOT NULL PRIMARY KEY,
    title   TEXT,
    author  TEXT,
    journal TEXT
);

CREATE TABLE IF NOT EXISTS gp_lit_ref (
    gp_accession             TEXT    NOT NULL REFERENCES genome_property(accession) ON DELETE CASCADE,
    literature_reference_pmid INTEGER NOT NULL REFERENCES literature_reference(pmid) ON DELETE CASCADE,
    list_order               INTEGER NOT NULL
);

CREATE TABLE IF NOT EXISTS gp_database_link (
    gp_accession TEXT NOT NULL REFERENCES genome_property(accession) ON DELETE CASCADE,
    db_id        TEXT NOT NULL,
    db_link      TEXT NOT NULL,
    other_params TEXT,
    comment      TEXT
);
"""

# ---------------------------------------------------------------------------
# Flatfile parser
#
# Mirrors the Perl parsing in GenomePropertiesIO.pm:
#   - parseFlatfile() splits the file on "//" record separators
#   - parseDESC() parses each record's header, references, DB refs
#   - parseSteps() parses step blocks (SN/ID/DN/RQ/EV/TG lines)
# ---------------------------------------------------------------------------

# Allowed reference tag transitions (mirrors $refTags in Perl)
_REF_TAGS = {
    'RC': {'RC', 'RN'},
    'RN': {'RM'},
    'RM': {'RT'},
    'RT': {'RT', 'RA'},
    'RA': {'RA', 'RL'},
    'RL': {'RL'},
}

# DR database patterns (mirrors the if/elsif chain in Perl parseDESC)
_DR_PATTERNS = [
    (re.compile(r'^DR  (KEGG);\s+(\S+);$'),              ('db_id', 'db_link')),
    (re.compile(r'^DR  (EcoCyc);\s+(\S+);$'),            ('db_id', 'db_link')),
    (re.compile(r'^DR  (MetaCyc);\s+(\S+);$'),           ('db_id', 'db_link')),
    (re.compile(r'^DR  (IUBMB);\s(\S+);\s(\S+);$'),      ('db_id', 'db_link', 'other_params')),
    (re.compile(r'^DR  (URL);\s+(\S+);$'),                ('db_id', 'db_link')),
    (re.compile(r'^DR  (Complex Portal);\s+(CPX-\S+);$'), ('db_id', 'db_link')),
    (re.compile(r'^DR  (PDBe);\s+(\S{4});$'),             ('db_id', 'db_link')),
]


def parse_flatfile(path):
    """Split the flatfile on '//' record separators and parse each record.

    This mirrors the Perl parseFlatfile() which sets $/ = "//" to read
    one record at a time, then feeds each chunk to parseDESC().
    """
    with open(path, 'r') as fh:
        content = fh.read()

    # Split on '//' exactly as Perl does with $/ = "//"
    records = content.split('//')

    properties = []
    for record in records:
        # Strip leading blank lines (mirrors: shift(@file) if $file[0] eq "")
        lines = record.split('\n')
        while lines and lines[0].strip() == '':
            lines.pop(0)
        # Skip empty records
        if not lines:
            continue
        prop = _parse_desc(lines)
        if prop and prop['AC'] is not None:
            properties.append(prop)

    return properties


def _parse_desc(lines):
    """Parse a single DESC record (list of lines without the '//' terminator).

    Mirrors parseDESC() in GenomePropertiesIO.pm.
    """
    params = {
        'AC': None, 'DE': None, 'TP': None, 'AU': None, 'TH': 0,
        'CC': None, 'private': None,
        'REFS': [], 'DBREFS': [], 'STEPS': [],
    }

    i = 0
    n = len(lines)
    while i < n:
        line = lines[i].rstrip('\n')

        # Single-value header tags: AC, DE, AU, TP, TH
        # Perl: $file[$i] =~ /^(AC|DE|AU|TP|TH|)\s{2}(.*)$/
        m = re.match(r'^(AC|DE|AU|TP|TH)\s{2}(.*)$', line)
        if m:
            tag, val = m.group(1), m.group(2)
            params[tag] = val
            i += 1
            continue

        # Private annotation: **  text
        m = re.match(r'^\*\*\s{2}(.*)$', line)
        if m:
            text = m.group(1)
            params['private'] = (params['private'] + ' ' + text) if params['private'] else text
            i += 1
            continue

        # Parent pointer: PN  GenPropNNNN (stored but not in DB schema)
        if re.match(r'^PN\s{2}GenProp\d{4}$', line):
            i += 1
            continue

        # Comment: CC  text
        m = re.match(r'^CC\s{2}(.*)$', line)
        if m:
            cc = m.group(1)
            params['CC'] = (params['CC'] + ' ' + cc) if params['CC'] else cc
            i += 1
            continue

        # Reference block: starts with RN or RC
        if re.match(r'^R[NC]\s{2}', line):
            ref, i = _parse_reference(lines, i)
            params['REFS'].append(ref)
            continue

        # Database reference block: starts with DC or DR
        if re.match(r'^D[CR]\s{2}', line):
            i = _parse_dbrefs(lines, i, params['DBREFS'])
            continue

        # Step separator '--' : start of steps section
        if line == '--':
            i += 1
            steps, i = _parse_steps(lines, i)
            params['STEPS'] = steps
            continue

        # Blank lines
        if line.strip() == '':
            i += 1
            continue

        # Anything else
        print(f"WARNING: unparsed line in {params.get('AC', '?')}: |{line}|",
              file=sys.stderr)
        i += 1

    return params


def _parse_reference(lines, i):
    """Parse an RN..RL reference block.

    Mirrors the REFLINE loop in Perl parseDESC(). Uses _REF_TAGS to know
    which tag transitions are valid, stopping when we hit a non-ref tag.
    """
    ref = {}
    n = len(lines)

    while i < n:
        line = lines[i].rstrip('\n')
        m = re.match(r'^(\w{2})\s{2}(.*)', line)
        if not m:
            break

        tag, val = m.group(1), m.group(2)
        # Only accept known reference tags
        if tag not in _REF_TAGS and tag not in ('RN', 'RC', 'RM', 'RT', 'RA', 'RL'):
            break

        if tag in ref:
            ref[tag] += ' ' + val
        else:
            ref[tag] = val
        i += 1

        # Check if the NEXT line is a valid continuation
        if i < n:
            next_line = lines[i].rstrip('\n')
            next_m = re.match(r'^(\S{2})', next_line)
            if next_m:
                next_tag = next_m.group(1)
                # Valid continuations from current tag
                if tag in _REF_TAGS and next_tag in _REF_TAGS.get(tag, set()):
                    continue
                # RL can be followed by another reference (RN/RC) or a non-ref tag
                if tag == 'RL':
                    break
                # If next_tag is a known ref tag, continue
                if next_tag in ('RN', 'RC', 'RM', 'RT', 'RA', 'RL'):
                    continue
            break

    # Clean reference number: strip brackets
    if 'RN' in ref:
        ref['RN'] = re.sub(r'[\[\]]', '', ref['RN'])

    return ref, i


def _parse_dbrefs(lines, i, dbrefs_list):
    """Parse a DC/DR block, appending entries to dbrefs_list.

    Mirrors the elsif($file[$i] =~ /^D\\w\\s{2}/) section in parseDESC().
    """
    n = len(lines)

    while i < n:
        line = lines[i].rstrip('\n')

        # Collect optional DC comment lines
        dc_comment = None
        while i < n:
            line = lines[i].rstrip('\n')
            m = re.match(r'^DC\s{2}(.*)', line)
            if m:
                text = m.group(1)
                dc_comment = (dc_comment + ' ' + text) if dc_comment else text
                i += 1
            else:
                break

        if i >= n:
            break

        line = lines[i].rstrip('\n')

        # Try each known DR pattern
        matched = False
        for pattern, fields in _DR_PATTERNS:
            m = pattern.match(line)
            if m:
                entry = {f: m.group(j + 1) for j, f in enumerate(fields)}
                if dc_comment:
                    entry['db_comment'] = dc_comment
                dbrefs_list.append(entry)
                matched = True
                i += 1
                break

        if not matched:
            # Unknown DR or non-DR line: check if it's a DR line we don't recognise
            if re.match(r'^DR', line):
                print(f"WARNING: unknown database reference: {line}", file=sys.stderr)
                i += 1
            else:
                # Not a DR line at all — back out (mirrors the else{$i--; last} in Perl)
                break

    return i


def _parse_steps(lines, i):
    """Parse step blocks until end of record.

    Mirrors parseSteps() in GenomePropertiesIO.pm.
    """
    n = len(lines)
    steps = []
    step = {}

    while i < n:
        line = lines[i].rstrip('\n')

        # Single-value step tags: SN, ID, DN, EC, RQ
        m = re.match(r'^(SN|ID|DN|EC|RQ)\s{2}(.*)$', line)
        if m:
            step[m.group(1)] = m.group(2)
            i += 1
            continue

        # IPR evidence with sufficiency: EV  IPR006219; TIGR00034; sufficient;
        m = re.match(r'^EV\s{2}(IPR\d{6});\s(\S+);\s(\S+);$', line)
        if m:
            ev = {'ipr': m.group(1), 'sig': m.group(2), 'sc': m.group(3), 'go': []}
            i += 1
            # Collect following TG lines
            while i < n:
                tg_m = re.match(r'^TG\s{2}(GO:\d+)', lines[i].rstrip('\n'))
                if tg_m:
                    ev['go'].append(tg_m.group(1))
                    i += 1
                else:
                    break
            step.setdefault('EVID', []).append(ev)
            continue

        # IPR evidence without sufficiency: EV  IPR002480; TIGR01358;
        m = re.match(r'^EV\s{2}(IPR\d{6});\s(\S+);$', line)
        if m:
            ev = {'ipr': m.group(1), 'sig': m.group(2), 'go': []}
            i += 1
            while i < n:
                tg_m = re.match(r'^TG\s{2}(GO:\d+)', lines[i].rstrip('\n'))
                if tg_m:
                    ev['go'].append(tg_m.group(1))
                    i += 1
                else:
                    break
            step.setdefault('EVID', []).append(ev)
            continue

        # GenProp evidence: EV  GenProp0791;
        m = re.match(r'^EV\s{2}(GenProp\d{4});$', line)
        if m:
            ev = {'gp': m.group(1), 'go': []}
            i += 1
            while i < n:
                tg_m = re.match(r'^TG\s{2}(GO:\d+)', lines[i].rstrip('\n'))
                if tg_m:
                    ev['go'].append(tg_m.group(1))
                    i += 1
                else:
                    break
            step.setdefault('EVID', []).append(ev)
            continue

        # Step separator '--'
        if line == '--':
            steps.append(step)
            step = {}
            i += 1
            continue

        # Blank line
        if line.strip() == '':
            i += 1
            continue

        # Anything else — warn and skip
        print(f"WARNING: unparsed step line: |{line}|", file=sys.stderr)
        i += 1

    # Last step (record ends without '//' since it was already split off)
    if step:
        steps.append(step)

    return steps, i


# ---------------------------------------------------------------------------
# GO term lookup (QuickGO REST API)
# ---------------------------------------------------------------------------

GO_API_URL = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{go_id}"

_GO_CACHE = {}


def fetch_go_term(go_id, session):
    """Return (term_text, category) for *go_id*."""
    if go_id in _GO_CACHE:
        return _GO_CACHE[go_id]

    url = GO_API_URL.format(go_id=go_id)
    for attempt in range(3):
        try:
            resp = session.get(url, headers={"Accept": "application/json"}, timeout=15)
            if resp.status_code == 200:
                data = resp.json()
                results = data.get("results", [])
                if results:
                    term = results[0].get("name", "")
                    aspect = results[0].get("aspect", "obsolete")
                    _GO_CACHE[go_id] = (term, aspect)
                    return term, aspect
                _GO_CACHE[go_id] = ("", "obsolete")
                return "", "obsolete"
            elif resp.status_code == 429:
                time.sleep(2 ** attempt)
                continue
            else:
                print(f"WARNING: GO lookup failed for {go_id} (HTTP {resp.status_code})",
                      file=sys.stderr)
                _GO_CACHE[go_id] = ("", "unknown")
                return "", "unknown"
        except requests.RequestException as exc:
            print(f"WARNING: GO lookup error for {go_id}: {exc}", file=sys.stderr)
            time.sleep(2 ** attempt)

    _GO_CACHE[go_id] = ("", "unknown")
    return "", "unknown"


# ---------------------------------------------------------------------------
# Database writer
#
# Mirrors gp2db(), addReferences2DB(), addDBrefs2DB(), addSteps2DB()
# in GenomePropertiesIO.pm
# ---------------------------------------------------------------------------

def create_database(db_path, properties, skip_go_lookup=False):
    """Create a SQLite database at *db_path*."""
    if os.path.exists(db_path):
        os.remove(db_path)

    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA journal_mode=WAL")
    cur = conn.cursor()
    cur.executescript(SCHEMA_SQL)
    # Disable FK checks during bulk loading — gp_step_evidence_gp may reference
    # a genome_property not yet inserted (the original Perl script used retries
    # for this). We verify FK integrity after all inserts are done.
    conn.execute("PRAGMA foreign_keys=OFF")

    http_session = requests.Session() if not skip_go_lookup else None

    # Collect all unique GO ids for pre-fetching
    go_ids_needed = set()
    for prop in properties:
        for step in prop['STEPS']:
            for ev in step.get('EVID', []):
                for go_id in ev.get('go', []):
                    go_ids_needed.add(go_id)

    # Populate go_terms table
    if go_ids_needed:
        if http_session:
            print(f"Fetching {len(go_ids_needed)} unique GO terms from QuickGO ...",
                  file=sys.stderr)
            for idx, go_id in enumerate(sorted(go_ids_needed), 1):
                term, category = fetch_go_term(go_id, http_session)
                cur.execute(
                    "INSERT OR IGNORE INTO go_terms (go_id, term, category) VALUES (?, ?, ?)",
                    (go_id, term, category),
                )
                if idx % 50 == 0:
                    print(f"  ... {idx}/{len(go_ids_needed)} GO terms fetched",
                          file=sys.stderr)
                    conn.commit()
            conn.commit()
            print(f"  ... done ({len(go_ids_needed)} GO terms).", file=sys.stderr)
        else:
            # Insert placeholder rows so FK constraints are satisfied
            print(f"Inserting {len(go_ids_needed)} GO term placeholders (--skip-go-lookup) ...",
                  file=sys.stderr)
            for go_id in sorted(go_ids_needed):
                cur.execute(
                    "INSERT OR IGNORE INTO go_terms (go_id, term, category) VALUES (?, ?, ?)",
                    (go_id, "", ""),
                )
            conn.commit()

    # Insert properties (mirrors gp2db loop)
    for prop in properties:
        acc = prop['AC']

        # genome_property row (mirrors update_or_create in gp2db)
        cur.execute(
            """INSERT OR REPLACE INTO genome_property
               (accession, description, type, author, threshold,
                comment, private, ispublic, checked)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            (acc, prop['DE'], prop['TP'], prop.get('AU', ''),
             int(prop.get('TH', 0)),
             prop.get('CC'),
             prop.get('private'),
             0,   # ispublic: not in flatfile, default 0
             1),  # checked: hardcoded to 1 as in Perl
        )

        # Literature references (mirrors addReferences2DB)
        for ref in prop['REFS']:
            pmid_str = ref.get('RM')
            if pmid_str is not None:
                pmid = int(pmid_str)
                cur.execute(
                    """INSERT OR REPLACE INTO literature_reference
                       (pmid, title, author, journal) VALUES (?, ?, ?, ?)""",
                    (pmid, ref.get('RT'), ref.get('RA'), ref.get('RL')),
                )
                rn = int(ref.get('RN', 0))
                cur.execute(
                    """INSERT INTO gp_lit_ref
                       (gp_accession, literature_reference_pmid, list_order)
                       VALUES (?, ?, ?)""",
                    (acc, pmid, rn),
                )

        # Database links (mirrors addDBrefs2DB)
        for dr in prop['DBREFS']:
            cur.execute(
                """INSERT INTO gp_database_link
                   (gp_accession, db_id, db_link, other_params, comment)
                   VALUES (?, ?, ?, ?, ?)""",
                (acc, dr['db_id'], dr['db_link'],
                 dr.get('other_params'), dr.get('db_comment')),
            )

        # Steps and evidence (mirrors addSteps2DB)
        for step in prop['STEPS']:
            sn = int(step.get('SN', 0))
            step_id = step.get('ID', '')
            display_name = step.get('DN')
            required = int(step.get('RQ', 0))

            cur.execute(
                """INSERT INTO gp_step
                   (gp_accession, step_number, step_id, step_display_name, required)
                   VALUES (?, ?, ?, ?, ?)""",
                (acc, sn, step_id, display_name, required),
            )
            auto_step = cur.lastrowid

            for ev in step.get('EVID', []):
                gos = ev.get('go', [])

                if 'gp' in ev:
                    # GenProp evidence (mirrors GpStepEvidenceGp create)
                    cur.execute(
                        """INSERT INTO gp_step_evidence_gp
                           (auto_step, gp_accession) VALUES (?, ?)""",
                        (auto_step, ev['gp']),
                    )
                    auto_gp_step = cur.lastrowid
                    for go_id in gos:
                        cur.execute(
                            "INSERT INTO gp_step_to_go (auto_gp_step, go_id) VALUES (?, ?)",
                            (auto_gp_step, go_id),
                        )
                else:
                    # InterPro evidence (mirrors GpStepEvidenceIpr create)
                    sufficient = 1 if ev.get('sc') == 'sufficient' else 0
                    cur.execute(
                        """INSERT INTO gp_step_evidence_ipr
                           (auto_step, interpro_acc, signature_acc, sufficient)
                           VALUES (?, ?, ?, ?)""",
                        (auto_step, ev['ipr'], ev['sig'], sufficient),
                    )
                    auto_ipr_step = cur.lastrowid
                    for go_id in gos:
                        cur.execute(
                            "INSERT INTO ipr_step_to_go (auto_ipr_step, go_id) VALUES (?, ?)",
                            (auto_ipr_step, go_id),
                        )

    conn.commit()

    # Verify FK integrity now that all data is loaded
    conn.execute("PRAGMA foreign_keys=ON")
    violations = conn.execute("PRAGMA foreign_key_check").fetchall()
    if violations:
        print(f"WARNING: {len(violations)} foreign key violations detected:",
              file=sys.stderr)
        for v in violations[:20]:
            print(f"  table={v[0]} rowid={v[1]} parent={v[2]} fkid={v[3]}",
                  file=sys.stderr)

    conn.close()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Parse a Genome Properties flatfile and create a SQLite database.",
    )
    parser.add_argument(
        "flatfile",
        help="Path to genomeProperties.txt",
    )
    parser.add_argument(
        "-o", "--output",
        default="genome_properties.db",
        help="Output SQLite database path (default: genome_properties.db)",
    )
    parser.add_argument(
        "--skip-go-lookup",
        action="store_true",
        help="Skip fetching GO term labels from the EBI API (go_terms table will be empty)",
    )
    args = parser.parse_args()

    print(f"Parsing {args.flatfile} ...", file=sys.stderr)
    properties = parse_flatfile(args.flatfile)
    print(f"Parsed {len(properties)} genome properties.", file=sys.stderr)

    print(f"Writing SQLite database to {args.output} ...", file=sys.stderr)
    create_database(args.output, properties, skip_go_lookup=args.skip_go_lookup)

    # Print summary
    conn = sqlite3.connect(args.output)
    cur = conn.cursor()
    tables = [
        'genome_property', 'gp_step', 'gp_step_evidence_ipr',
        'gp_step_evidence_gp', 'go_terms', 'gp_step_to_go',
        'ipr_step_to_go', 'literature_reference', 'gp_lit_ref',
        'gp_database_link',
    ]
    print("\nDatabase summary:", file=sys.stderr)
    for table in tables:
        count = cur.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
        print(f"  {table:30s} {count:>6d} rows", file=sys.stderr)
    conn.close()
    print(f"\nDone. Database written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
