#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import asyncio
import aiohttp
import csv
import re
import json
from tqdm import tqdm

ALPHAFOLD_API = "https://www.alphafold.ebi.ac.uk/api/prediction"
MGNIFY_WEB = "https://www.ebi.ac.uk/metagenomics/proteins"
RCSB_PDB_API = "https://data.rcsb.org/graphql"

CONCURRENT_REQUESTS = 20


async def fetch_alphafold_description(session, uniprot_id, semaphore, retries=3):
    url = ALPHAFOLD_API + "/" + uniprot_id

    for attempt in range(retries):
        try:
            async with semaphore:
                async with session.get(url, timeout=30) as r:
                    if r.status != 200:
                        return uniprot_id, "HTTP {}".format(r.status)

                    data = await r.json()
                    if isinstance(data, list) and data:
                        return uniprot_id, data[0].get("uniprotDescription", "No description found")
                    return uniprot_id, "No data"
        except Exception as e:
            if attempt == retries - 1:
                return uniprot_id, "Failed: {}".format(e)
            await asyncio.sleep(1)


async def fetch_pdb_annotations(session, pdb_id, semaphore, retries=3):
    query = {
        "query": """
        query($id: String!) {
          entry(entry_id: $id) {
            struct { title }
            polymer_entities {
              rcsb_polymer_entity { pdbx_description }
              pfams {
                rcsb_pfam_accession
                rcsb_pfam_identifier
                rcsb_pfam_description
              }
            }
          }
        }
        """,
        "variables": {"id": pdb_id.upper()}
    }

    for attempt in range(retries):
        try:
            async with semaphore:
                async with session.post(RCSB_PDB_API, json=query, timeout=30) as r:
                    if r.status == 404:
                        return pdb_id, {
                            "pdb_id": pdb_id,
                            "title": "Not found",
                            "pfam": [],
                            "description": "Not found"
                        }

                    if r.status != 200:
                        raise RuntimeError("HTTP {}".format(r.status))

                    raw = await r.json()
                    entry = raw.get("data", {}).get("entry")

                    ann = {
                        "pdb_id": pdb_id,
                        "title": "No title",
                        "pfam": [],
                        "description": "No description"
                    }

                    if not entry:
                        return pdb_id, ann

                    if entry.get("struct", {}).get("title"):
                        ann["title"] = entry["struct"]["title"]

                    for ent in entry.get("polymer_entities", []):
                        desc = ent.get("rcsb_polymer_entity", {}).get("pdbx_description")
                        if desc:
                            ann["description"] = desc

                        for pf in ent.get("pfams") or []:
                            ann["pfam"].append({
                                "accession": pf.get("rcsb_pfam_accession", ""),
                                "name": pf.get("rcsb_pfam_identifier", ""),
                                "description": pf.get("rcsb_pfam_description", "")
                            })

                    return pdb_id, ann

        except Exception as e:
            if attempt == retries - 1:
                return pdb_id, {
                    "pdb_id": pdb_id,
                    "title": "Failed",
                    "pfam": [],
                    "description": str(e),
                    "error": str(e)
                }
            await asyncio.sleep(2)


async def fetch_mgnify_from_web(session, mgyp_id, semaphore, retries=3):
    url = MGNIFY_WEB + "/" + mgyp_id + "/"

    for attempt in range(retries):
        try:
            async with semaphore:
                async with session.get(url, timeout=30) as r:
                    if r.status == 404:
                        return mgyp_id, {
                            "mgyp_id": mgyp_id,
                            "sequence_length": "N/A",
                            "pfam": [],
                            "interpro": [],
                            "go_terms": [],
                            "description": "Not found"
                        }

                    if r.status != 200:
                        raise RuntimeError("HTTP {}".format(r.status))

                    html = await r.text()

                    ann = {
                        "mgyp_id": mgyp_id,
                        "sequence_length": "N/A",
                        "pfam": [],
                        "interpro": [],
                        "go_terms": [],
                        "description": "No description"
                    }

                    seq = re.search(r'<div id="proteinSequenceContainer".*?>(.*?)</div>', html, re.S)
                    if seq:
                        ann["sequence_length"] = len(re.sub(r"\s+", "", seq.group(1)))

                    pfam_json = re.search(
                        r'<script id="pfam-annotations-data".*?>(\[.*?\])</script>',
                        html,
                        re.S
                    )

                    if pfam_json:
                        for p in json.loads(pfam_json.group(1)):
                            ann["pfam"].append({
                                "accession": p.get("accession", ""),
                                "name": p.get("name", ""),
                                "description": ""
                            })

                    if not ann["pfam"]:
                        rows = re.findall(r'(PF\d{5}).*?<td.*?>([^<]+)</td>', html, re.S)
                        for pf, name in rows:
                            ann["pfam"].append({
                                "accession": pf,
                                "name": name.strip(),
                                "description": ""
                            })

                    iprs = re.findall(r'(IPR\d{6}).*?<td.*?>([^<]+)</td>', html)
                    for ipr, name in iprs:
                        ann["interpro"].append({
                            "accession": ipr,
                            "name": name.strip(),
                            "type": ""
                        })

                    gos = re.findall(r'(GO:\d{7}).*?<td.*?>([^<]+)</td>', html)
                    for go, name in gos:
                        ann["go_terms"].append({
                            "accession": go,
                            "name": name.strip(),
                            "category": ""
                        })

                    return mgyp_id, ann

        except Exception as e:
            if attempt == retries - 1:
                return mgyp_id, {
                    "mgyp_id": mgyp_id,
                    "sequence_length": "N/A",
                    "pfam": [],
                    "interpro": [],
                    "go_terms": [],
                    "description": str(e),
                    "error": str(e)
                }
            await asyncio.sleep(2)


def extract_unique_ids(m8_file, database):
    ids = set()

    with open(m8_file) as f:
        for line in f:
            cols = line.rstrip().split("\t")
            if len(cols) < 2:
                continue

            target = cols[1]

            if database == "alphafold" and target.startswith("AF-") and "-F1" in target:
                ids.add(target.split("-F1")[0].replace("AF-", ""))

            elif database == "mgnify" and "MGYP" in target:
                mg = target[target.find("MGYP"):]
                mg = mg.split(".")[0].split("_")[0]
                ids.add(mg)

            elif database == "pdb":
                pid = target.split("-")[0].split(".")[0].lower()
                if len(pid) == 4:
                    ids.add(pid)

    return sorted(ids)


async def fetch_all_annotations(ids, database):
    sem = asyncio.Semaphore(CONCURRENT_REQUESTS)
    conn = aiohttp.TCPConnector(limit=CONCURRENT_REQUESTS)

    async with aiohttp.ClientSession(connector=conn) as session:
        if database == "alphafold":
            tasks = [fetch_alphafold_description(session, i, sem) for i in ids]
        elif database == "mgnify":
            tasks = [fetch_mgnify_from_web(session, i, sem) for i in ids]
        else:
            tasks = [fetch_pdb_annotations(session, i, sem) for i in ids]

        results = {}
        for f in tqdm(asyncio.as_completed(tasks), total=len(tasks)):
            k, v = await f
            results[k] = v

        return results


def format_ann(lst):
    if not lst:
        return "None"
    return "; ".join("{} ({})".format(x["accession"], x["name"]) for x in lst)


def merge_alphafold_annotations(m8, ann, out):
    with open(m8) as inp, open(out, "w", newline="") as outp:
        w = csv.writer(outp, delimiter="\t")
        w.writerow([
            "Query_ID", "Target_ID", "UniProt_ID", "Identity", "Length",
            "Mismatches", "GapOpen", "Q_start", "Q_end", "S_start", "S_end",
            "E-value", "BitScore", "Description"
        ])

        for line in inp:
            p = line.rstrip().split("\t")
            if len(p) < 12:
                continue

            tid = p[1]
            uid = tid.split("-F1")[0].replace("AF-", "")
            desc = ann.get(uid, "No description found")
            w.writerow(p[:1] + [tid, uid] + p[2:] + [desc])


def merge_mgnify_annotations(m8, ann, out):
    with open(m8) as inp, open(out, "w", newline="") as outp:
        w = csv.writer(outp, delimiter="\t")
        w.writerow([
            "Query_ID", "Target_ID", "MGYP_ID", "Identity", "Length",
            "Mismatches", "GapOpen", "Q_start", "Q_end", "S_start", "S_end",
            "E-value", "BitScore", "Seq_Length", "Pfam_Annotations",
            "InterPro_Annotations", "GO_Terms", "Description"
        ])

        for line in inp:
            p = line.rstrip().split("\t")
            if len(p) < 12:
                continue

            t = p[1]
            mg = t[t.find("MGYP"):] if "MGYP" in t else t
            mg = mg.split(".")[0].split("_")[0]

            a = ann.get(mg, {})
            w.writerow(
                p[:1] + [t, mg] + p[2:] + [
                    a.get("sequence_length", "N/A"),
                    format_ann(a.get("pfam", [])),
                    format_ann(a.get("interpro", [])),
                    format_ann(a.get("go_terms", [])),
                    a.get("description", "No description")
                ]
            )


def merge_pdb_annotations(m8, ann, out):
    with open(m8) as inp, open(out, "w", newline="") as outp:
        w = csv.writer(outp, delimiter="\t")
        w.writerow([
            "Query_ID", "Target_ID", "PDB_ID", "Identity", "Length",
            "Mismatches", "GapOpen", "Q_start", "Q_end", "S_start", "S_end",
            "E-value", "BitScore", "Pfam_Annotations", "Title", "Description"
        ])

        for line in inp:
            p = line.rstrip().split("\t")
            if len(p) < 12:
                continue

            t = p[1]
            pid = t.split("-")[0].split(".")[0].lower()
            a = ann.get(pid, {})
            w.writerow(
                p[:1] + [t, pid] + p[2:] + [
                    format_ann(a.get("pfam", [])),
                    a.get("title", "No title"),
                    a.get("description", "No description")
                ]
            )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("-d", "--database", choices=["alphafold", "mgnify", "pdb"], required=True)
    ap.add_argument("-c", "--concurrent", type=int, default=20)
    args = ap.parse_args()

    global CONCURRENT_REQUESTS
    CONCURRENT_REQUESTS = args.concurrent

    ids = extract_unique_ids(args.input, args.database)
    if not ids:
        raise SystemExit("No valid IDs found")

    loop = asyncio.get_event_loop()
    ann = loop.run_until_complete(fetch_all_annotations(ids, args.database))

    if args.database == "alphafold":
        merge_alphafold_annotations(args.input, ann, args.output)
    elif args.database == "mgnify":
        merge_mgnify_annotations(args.input, ann, args.output)
    else:
        merge_pdb_annotations(args.input, ann, args.output)


if __name__ == "__main__":
    main()
