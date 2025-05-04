import argparse
import asyncio
import aiohttp
import csv
import os
from tqdm import tqdm

API_BASE = "https://www.alphafold.ebi.ac.uk/api/prediction"
CONCURRENT_REQUESTS = 50  # You can adjust this based on your internet and server limit
async def fetch_description(session, uniprot_id, semaphore, retries=3):
    url = f"{API_BASE}/{uniprot_id}"
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
                        return uniprot_id, f"Error {response.status}"
        except Exception as e:
            if attempt == retries - 1:
                return uniprot_id, f"Failed: {e}"
            await asyncio.sleep(1)

def extract_unique_uniprot_ids(m8_file):
    uniprot_ids = set()
    with open(m8_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                target = parts[1]
                # Extract UniProt ID from AF-XXXX-F1-model_vX
                if target.startswith("AF-") and "-F1" in target:
                    id_part = target.split("-F1")[0]
                    uniprot_id = id_part.replace("AF-", "")
                    uniprot_ids.add(uniprot_id)
    return sorted(uniprot_ids)

async def fetch_all_descriptions(uniprot_ids):
    semaphore = asyncio.Semaphore(CONCURRENT_REQUESTS)
    connector = aiohttp.TCPConnector(limit=CONCURRENT_REQUESTS)
    timeout = aiohttp.ClientTimeout(total=60)
    async with aiohttp.ClientSession(connector=connector, timeout=timeout) as session:
        tasks = []
        for uniprot_id in uniprot_ids:
            tasks.append(fetch_description(session, uniprot_id, semaphore))
        results = []
        for f in tqdm(asyncio.as_completed(tasks), total=len(tasks), desc="Fetching annotations"):
            result = await f
            results.append(result)
        return dict(results)  # Return a dictionary for easy lookup

def merge_annotations_with_m8(m8_file, descriptions, output_file):
    with open(m8_file) as infile, open(output_file, "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow(["Query_ID", "UniProt_ID", "Identity", "Length", "Mismatches", "GapOpen", "Q_start", "Q_end", "S_start", "S_end", "E-value", "BitScore", "TaxID", "Description"])

        for line in infile:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                query_id = parts[0]
                uniprot_id = parts[1].split("-F1")[0].replace("AF-", "")
                description = descriptions.get(uniprot_id, "No description found")
                writer.writerow(parts + [description])

def main():
    parser = argparse.ArgumentParser(description="Annotate Foldseek AlphaFold matches using AlphaFold API")
    parser.add_argument("-i", "--input", required=True, help="Input Foldseek .m8 file")
    parser.add_argument("-o", "--output", required=True, help="Output annotation TSV file")
    args = parser.parse_args()

    print("Extracting UniProt IDs from m8 file...")
    uniprot_ids = extract_unique_uniprot_ids(args.input)
    print(f"Found {len(uniprot_ids)} unique UniProt IDs.")

    try:
        # Fetch descriptions for the UniProt IDs
        descriptions = asyncio.run(fetch_all_descriptions(uniprot_ids))
        
        # Merge the annotations with the original .m8 file and write to the output file
        merge_annotations_with_m8(args.input, descriptions, args.output)
        
        print(f"\n? Annotation completed. Output written to: {args.output}")
    except Exception as e:
        print(f"\n? Error during annotation: {e}")

if __name__ == "__main__":
    main()
