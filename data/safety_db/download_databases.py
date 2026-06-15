#!/usr/bin/env python3
"""
Download safety screening databases for TOPE_DEEP N6 local implementation.

Run at Docker build time:
  python data/safety_db/download_databases.py

Databases downloaded:
  1. AllergenOnline v19 (University of Nebraska-Lincoln)
     Source: allergenonline.org - WHO regulatory allergen database
     Size: ~2MB
     License: free for non-commercial use

  2. Human UniProt Swiss-Prot reviewed (taxon 9606)
     Source: uniprot.org - experimentally validated human proteins
     Size: ~45MB (FASTA)
     License: CC BY 4.0

Both are cached in data/safety_db/ and baked into the Docker image.
Update quarterly by re-running this script and rebuilding.
"""

import os
import sys
import requests
import datetime

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)))
os.makedirs(OUT_DIR, exist_ok=True)

VERSION = datetime.datetime.now().strftime("%Y.Q") + str(
    (datetime.datetime.now().month - 1) // 3 + 1
)


def download(url: str, dest: str, description: str) -> bool:
    print(f"Downloading {description}...")
    print(f"  URL: {url}")
    print(f"  Destination: {dest}")
    try:
        resp = requests.get(url, stream=True, timeout=120)
        resp.raise_for_status()
        total = int(resp.headers.get("content-length", 0))
        downloaded = 0
        with open(dest, "wb") as f:
            for chunk in resp.iter_content(chunk_size=8192):
                f.write(chunk)
                downloaded += len(chunk)
                if total:
                    pct = downloaded / total * 100
                    print(f"\r  {pct:.1f}%", end="", flush=True)
        print(f"\n  Done: {downloaded / 1024 / 1024:.1f} MB")
        return True
    except Exception as e:
        print(f"\n  FAILED: {e}")
        return False


def main():
    success = True

    # 1. AllergenOnline - FASTA of all allergen sequences
    allergen_dest = os.path.join(OUT_DIR, "allergenonline_allergens.fasta")
    if not os.path.exists(allergen_dest):
        ok = download(
            url="https://www.allergenonline.org/download/fasta_download.php",
            dest=allergen_dest,
            description="AllergenOnline allergen sequences (WHO regulatory database)",
        )
        if not ok:
            # Fallback: direct FASTA download URL
            ok = download(
                url="https://www.allergenonline.org/download/allergenonline_allergens_uniprot.fasta",
                dest=allergen_dest,
                description="AllergenOnline (fallback URL)",
            )
        success = success and ok
    else:
        print(f"AllergenOnline already exists: {allergen_dest}")

    # 2. Human Swiss-Prot reviewed
    human_dest = os.path.join(OUT_DIR, "human_swissprot_reviewed.fasta")
    if not os.path.exists(human_dest):
        ok = download(
            url="https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28reviewed%3Atrue%29+AND+%28organism_id%3A9606%29",
            dest=human_dest,
            description="Human UniProt Swiss-Prot reviewed (taxon 9606)",
        )
        success = success and ok
    else:
        print(f"Human Swiss-Prot already exists: {human_dest}")

    # Write version file
    version_file = os.path.join(OUT_DIR, "VERSION")
    with open(version_file, "w") as f:
        f.write(f"ALLERGENONLINE_VERSION={VERSION}\n")
        f.write(f"HUMAN_SWISSPROT_VERSION={VERSION}\n")
        f.write(f"DOWNLOADED={datetime.datetime.now().isoformat()}\n")
    print(f"\nVersion recorded: {VERSION}")

    if success:
        print("\nAll databases downloaded successfully.")
        print("Set environment variables in Railway:")
        print(f"  ALLERGENONLINE_VERSION={VERSION}")
        print(f"  HUMAN_SWISSPROT_VERSION={VERSION}")
    else:
        print("\nSome downloads failed. Check URLs above.")
        sys.exit(1)


if __name__ == "__main__":
    main()