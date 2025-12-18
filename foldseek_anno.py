import argparse
import asyncio
import aiohttp
import csv
import re
import json
from tqdm import tqdm

# API Configuration
ALPHAFOLD_API = "https://www.alphafold.ebi.ac.uk/api/prediction"
MGNIFY_WEB = "https://www.ebi.ac.uk/metagenomics/proteins"
CONCURRENT_REQUESTS = 20
VERBOSE = False

async def fetch_alphafold_description(session, uniprot_id, semaphore, retries=3):
    """Fetch protein description from AlphaFold API"""
    url = "{}/{}".format(ALPHAFOLD_API, uniprot_id)
    for attempt in range(retries):
        try:
            async with semaphore:
                async with session.get(url, timeout=30) as response:
                    if response.status == 200:
                        data = await response.json()
                        if isinstance(data, list) and data:
                            return uniprot_id, data[0].get("uniprotDescription", "No description found")
                        else:
                            return uniprot_id, "No data"
                    else:
                        return uniprot_id, "Error {}".format(response.status)
        except Exception as e:
            if attempt == retries - 1:
                return uniprot_id, "Failed: {}".format(e)
            await asyncio.sleep(1)

async def fetch_mgnify_from_web(session, mgyp_id, semaphore, retries=3):
    """Fetch functional annotations by scraping MGnify website"""
    url = "{}/{}/".format(MGNIFY_WEB, mgyp_id)
    
    if VERBOSE:
        print("Fetching: {}".format(url))
    
    for attempt in range(retries):
        try:
            async with semaphore:
                async with session.get(url, timeout=30) as response:
                    status = response.status
                    
                    if VERBOSE:
                        print("  Status: {}".format(status))
                    
                    if status == 200:
                        html = await response.text()
                        
                        # Initialize annotation structure
                        annotations = {
                            'mgyp_id': mgyp_id,
                            'sequence_length': 'N/A',
                            'pfam': [],
                            'interpro': [],
                            'go_terms': [],
                            'description': 'No description'
                        }
                        
                        # Extract sequence length
                        seq_pattern = r'<div id="proteinSequenceContainer"[^>]*>([\s\S]*?)</div>'
                        seq_match = re.search(seq_pattern, html)
                        if seq_match:
                            sequence = seq_match.group(1).strip()
                            # Remove all whitespace from sequence
                            sequence = re.sub(r'\s+', '', sequence)
                            annotations['sequence_length'] = len(sequence)
                        
                        # Extract Pfam annotations
                        pfam_json_pattern = r'<script id="pfam-annotations-data" type="application/json">(\[.*?\])</script>'
                        pfam_json_match = re.search(pfam_json_pattern, html, re.DOTALL)
                        
                        if pfam_json_match:
                            try:
                                pfam_data = json.loads(pfam_json_match.group(1))
                                for pfam_entry in pfam_data:
                                    annotations['pfam'].append({
                                        'accession': pfam_entry.get('accession', ''),
                                        'name': pfam_entry.get('name', ''),
                                        'description': ''
                                    })
                            except:
                                pass
                        
                        # If JSON parsing failed, fall back to table parsing
                        if not annotations['pfam']:
                            # Parse Pfam table rows
                            pfam_row_pattern = r'<tr class="vf-table__row pfam-item"[^>]*>.*?<a[^>]*>(PF\d{5})</a>.*?<td class="vf-table__cell">([^<]+)</td>'
                            pfam_rows = re.findall(pfam_row_pattern, html, re.DOTALL)
                            for pfam_id, pfam_name in pfam_rows:
                                annotations['pfam'].append({
                                    'accession': pfam_id,
                                    'name': pfam_name.strip(),
                                    'description': ''
                                })
                        
                        # Extract InterPro annotations from links
                        # Look for InterPro links in the page
                        ipr_pattern = r'interpro/entry/[^/]+/(IPR\d{6})"[^>]*>IPR\d{6}</a>\s*</td>\s*<td[^>]*>([^<]+)</td>'
                        ipr_matches = re.findall(ipr_pattern, html)
                        for ipr_id, ipr_name in ipr_matches:
                            annotations['interpro'].append({
                                'accession': ipr_id,
                                'name': ipr_name.strip(),
                                'type': ''
                            })
                        
                        # Extract GO terms
                        go_pattern = r'<a[^>]*>(GO:\d{7})</a>\s*</td>\s*<td[^>]*>([^<]+)</td>'
                        go_matches = re.findall(go_pattern, html)
                        for go_id, go_name in go_matches:
                            annotations['go_terms'].append({
                                'accession': go_id,
                                'name': go_name.strip(),
                                'category': ''
                            })
                        
                        if VERBOSE:
                            print("  Length: {}".format(annotations['sequence_length']))
                            print("  Pfam: {}, InterPro: {}, GO: {}".format(
                                len(annotations['pfam']),
                                len(annotations['interpro']),
                                len(annotations['go_terms'])
                            ))
                        
                        return mgyp_id, annotations
                    
                    elif status == 404:
                        return mgyp_id, {
                            'mgyp_id': mgyp_id,
                            'sequence_length': 'N/A',
                            'pfam': [],
                            'interpro': [],
                            'go_terms': [],
                            'description': 'Not found in MGnify'
                        }
                    else:
                        if attempt == retries - 1:
                            return mgyp_id, {
                                'mgyp_id': mgyp_id,
                                'sequence_length': 'N/A',
                                'pfam': [],
                                'interpro': [],
                                'go_terms': [],
                                'description': 'HTTP Error {}'.format(status),
                                'error': 'HTTP {}'.format(status)
                            }
        
        except Exception as e:
            if VERBOSE:
                print("  Exception: {}".format(e))
            if attempt == retries - 1:
                return mgyp_id, {
                    'mgyp_id': mgyp_id,
                    'sequence_length': 'N/A',
                    'pfam': [],
                    'interpro': [],
                    'go_terms': [],
                    'description': 'Failed: {}'.format(str(e)),
                    'error': str(e)
                }
            await asyncio.sleep(2)

def extract_unique_ids(m8_file, database):
    """Extract unique protein IDs from m8 file based on database type"""
    ids = set()
    
    with open(m8_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                target = parts[1]
                
                if database == 'alphafold':
                    if target.startswith("AF-") and "-F1" in target:
                        id_part = target.split("-F1")[0]
                        protein_id = id_part.replace("AF-", "")
                        ids.add(protein_id)
                
                elif database == 'mgnify':
                    if target.startswith("MGYP") or "MGYP" in target:
                        if target.startswith("MGYP"):
                            protein_id = target
                        else:
                            mgyp_start = target.find("MGYP")
                            protein_id = target[mgyp_start:]
                        
                        protein_id = protein_id.split('.')[0].split('_')[0].split()[0]
                        
                        if protein_id.startswith("MGYP") and len(protein_id) > 4:
                            ids.add(protein_id)
    
    return sorted(ids)

async def fetch_all_annotations(protein_ids, database):
    """Fetch annotations for all protein IDs from specified database"""
    semaphore = asyncio.Semaphore(CONCURRENT_REQUESTS)
    connector = aiohttp.TCPConnector(limit=CONCURRENT_REQUESTS)
    timeout = aiohttp.ClientTimeout(total=60)
    
    async with aiohttp.ClientSession(connector=connector, timeout=timeout) as session:
        tasks = []
        
        if database == 'alphafold':
            for protein_id in protein_ids:
                tasks.append(fetch_alphafold_description(session, protein_id, semaphore))
        elif database == 'mgnify':
            for protein_id in protein_ids:
                tasks.append(fetch_mgnify_from_web(session, protein_id, semaphore))
        
        results = []
        for f in tqdm(asyncio.as_completed(tasks), total=len(tasks), desc="Fetching {} annotations".format(database)):
            result = await f
            results.append(result)
        
        return dict(results)

def format_annotation_list(annotations, annotation_type='pfam'):
    """Format annotation list into readable string"""
    if not annotations:
        return "None"
    
    if annotation_type == 'pfam':
        return "; ".join(["{} ({})".format(a['accession'], a['name']) for a in annotations])
    elif annotation_type == 'interpro':
        return "; ".join(["{} ({})".format(a['accession'], a['name']) for a in annotations])
    elif annotation_type == 'go':
        return "; ".join(["{} ({})".format(a['accession'], a['name']) for a in annotations])
    
    return "None"

def merge_alphafold_annotations(m8_file, descriptions, output_file):
    """Merge AlphaFold annotations with m8 file"""
    with open(m8_file) as infile, open(output_file, "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow(["Query_ID", "Target_ID", "UniProt_ID", "Identity", "Length", 
                        "Mismatches", "GapOpen", "Q_start", "Q_end", "S_start", "S_end", 
                        "E-value", "BitScore", "Description"])
        
        for line in infile:
            parts = line.strip().split('\t')
            if len(parts) >= 12:
                target_id = parts[1]
                uniprot_id = target_id.split("-F1")[0].replace("AF-", "")
                description = descriptions.get(uniprot_id, "No description found")
                writer.writerow(parts[:1] + [target_id, uniprot_id] + parts[2:] + [description])

def merge_mgnify_annotations(m8_file, annotations, output_file):
    """Merge MGnify annotations with m8 file"""
    with open(m8_file) as infile, open(output_file, "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow(["Query_ID", "Target_ID", "MGYP_ID", "Identity", "Length", 
                        "Mismatches", "GapOpen", "Q_start", "Q_end", "S_start", "S_end", 
                        "E-value", "BitScore", "Seq_Length", "Pfam_Annotations", 
                        "InterPro_Annotations", "GO_Terms", "Description"])
        
        for line in infile:
            parts = line.strip().split('\t')
            if len(parts) >= 12:
                target_id = parts[1]
                
                if target_id.startswith("MGYP") or "MGYP" in target_id:
                    if target_id.startswith("MGYP"):
                        mgyp_id = target_id
                    else:
                        mgyp_start = target_id.find("MGYP")
                        mgyp_id = target_id[mgyp_start:]
                    
                    mgyp_id = mgyp_id.split('.')[0].split('_')[0].split()[0]
                else:
                    mgyp_id = target_id
                
                annot = annotations.get(mgyp_id, {})
                
                seq_length = annot.get('sequence_length', 'N/A')
                pfam_str = format_annotation_list(annot.get('pfam', []), 'pfam')
                interpro_str = format_annotation_list(annot.get('interpro', []), 'interpro')
                go_str = format_annotation_list(annot.get('go_terms', []), 'go')
                description = annot.get('description', 'No description')
                
                writer.writerow(parts[:1] + [target_id, mgyp_id] + parts[2:] + 
                              [seq_length, pfam_str, interpro_str, go_str, description])

def test_single_id(mgyp_id):
    """Test fetching a single MGnify ID"""
    print("\n=== Testing single ID: {} ===".format(mgyp_id))
    
    async def test():
        semaphore = asyncio.Semaphore(1)
        connector = aiohttp.TCPConnector(limit=1)
        timeout = aiohttp.ClientTimeout(total=60)
        
        async with aiohttp.ClientSession(connector=connector, timeout=timeout) as session:
            result = await fetch_mgnify_from_web(session, mgyp_id, semaphore)
            return result
    
    loop = asyncio.get_event_loop()
    result_id, result_data = loop.run_until_complete(test())
    
    print("\nResult for {}:".format(result_id))
    print("  Sequence Length: {}".format(result_data.get('sequence_length', 'N/A')))
    print("  Description: {}".format(result_data.get('description', 'N/A')))
    print("  Pfam annotations: {}".format(len(result_data.get('pfam', []))))
    if result_data.get('pfam'):
        for pfam in result_data['pfam']:
            print("    - {} ({})".format(pfam['accession'], pfam['name']))
    print("  InterPro annotations: {}".format(len(result_data.get('interpro', []))))
    if result_data.get('interpro'):
        for ipr in result_data['interpro'][:3]:  # Show first 3
            print("    - {} ({})".format(ipr['accession'], ipr['name']))
    print("  GO terms: {}".format(len(result_data.get('go_terms', []))))
    if 'error' in result_data:
        print("  ERROR: {}".format(result_data['error']))

def main():
    parser = argparse.ArgumentParser(
        description="Annotate Foldseek matches using AlphaFold API or MGnify website",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Annotate AlphaFold matches
  python %(prog)s -i foldseek_alphafold.m8 -o annotated_alphafold.tsv -d alphafold
  
  # Annotate MGnify matches by scraping website
  python %(prog)s -i foldseek_mgnify.m8 -o annotated_mgnify.tsv -d mgnify
  
  # Test a single MGnify ID with verbose output
  python %(prog)s -d mgnify --test MGYP000612992776 -v
        """
    )
    
    parser.add_argument("-i", "--input", 
                       help="Input Foldseek .m8 file")
    parser.add_argument("-o", "--output", 
                       help="Output annotation TSV file")
    parser.add_argument("-d", "--database", required=True, 
                       choices=['alphafold', 'mgnify'],
                       help="Database to query (alphafold or mgnify)")
    parser.add_argument("-c", "--concurrent", type=int, default=20,
                       help="Number of concurrent requests (default: 20)")
    parser.add_argument("-v", "--verbose", action='store_true',
                       help="Enable verbose output for debugging")
    parser.add_argument("--test", type=str,
                       help="Test mode: fetch annotations for a single protein ID")
    
    args = parser.parse_args()
    
    global CONCURRENT_REQUESTS, VERBOSE
    CONCURRENT_REQUESTS = args.concurrent
    VERBOSE = args.verbose
    
    if args.test:
        test_single_id(args.test)
        return
    
    if not args.input or not args.output:
        parser.error("--input and --output are required unless using --test mode")
    
    print("Extracting protein IDs from m8 file for {}...".format(args.database))
    protein_ids = extract_unique_ids(args.input, args.database)
    print("Found {} unique protein IDs.".format(len(protein_ids)))
    
    if VERBOSE and len(protein_ids) > 0:
        print("First few IDs: {}".format(protein_ids[:5]))
    
    if len(protein_ids) == 0:
        print("WARNING: No {} IDs found in the input file.".format(args.database))
        return
    
    try:
        loop = asyncio.get_event_loop()
        annotations = loop.run_until_complete(
            fetch_all_annotations(protein_ids, args.database)
        )
        
        if args.database == 'alphafold':
            merge_alphafold_annotations(args.input, annotations, args.output)
        elif args.database == 'mgnify':
            merge_mgnify_annotations(args.input, annotations, args.output)
        
        print("\nSUCCESS: Annotation completed. Output written to: {}".format(args.output))
        
        if args.database == 'mgnify':
            pfam_count = sum(1 for a in annotations.values() if a.get('pfam'))
            failed_count = sum(1 for a in annotations.values() if 'error' in a)
            print("\nSummary:")
            print("  Total proteins queried: {}".format(len(annotations)))
            print("  Proteins with Pfam annotations: {}".format(pfam_count))
            if failed_count > 0:
                print("  Proteins with errors: {}".format(failed_count))
        
    except Exception as e:
        print("\nERROR during annotation: {}".format(e))
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
