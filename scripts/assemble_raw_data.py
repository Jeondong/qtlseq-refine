"""Restore and verify the complete raw Data S3/S4 from repository-local parts."""
from pathlib import Path
import argparse,hashlib,json,lzma

ROOT=Path(__file__).resolve().parents[1]

def assemble_archives(root,output_dir,extract=False):
    metadata=json.loads((root/'provenance/raw_data_archives.json').read_text(encoding='utf-8'))
    output_dir.mkdir(parents=True,exist_ok=True);records=0
    for entry in metadata:
        target=output_dir/entry['name'];whole=hashlib.sha256();written=0
        with target.open('wb') as destination:
            for part in entry['parts']:
                data=(root/part['path']).read_bytes()
                if len(data)!=part['size'] or hashlib.sha256(data).hexdigest()!=part['sha256']:
                    raise ValueError(f'Part missing or changed: {part["path"]}')
                destination.write(data);whole.update(data);written+=len(data)
        if written!=entry['size'] or whole.hexdigest()!=entry['sha256']:
            raise ValueError(f'Reconstructed archive mismatch: {entry["name"]}')
        digest=hashlib.sha256();lines=0;size=0
        tsv=target.with_suffix('');destination=tsv.open('wb') if extract else None
        try:
            with lzma.open(target,'rb') as source:
                while block:=source.read(4*1024*1024):
                    digest.update(block);lines+=block.count(b'\n');size+=len(block)
                    if destination:destination.write(block)
        finally:
            if destination:destination.close()
        if digest.hexdigest()!=entry['uncompressed_sha256'] or size!=entry['uncompressed_bytes'] or lines!=entry['records']+1:
            raise ValueError(f'Decompressed TSV mismatch: {entry["name"]}')
        records+=entry['records']
        print(f'Data S{entry["number"]}: {entry["records"]:,} records verified; {tsv if extract else target}',flush=True)
    return dict(raw_archives_verified=len(metadata),raw_snp_records_verified=records,raw_tsv_hashes_match_original_submission=True)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=ROOT)
    parser.add_argument('--output-dir',type=Path,default=Path('outputs/raw_supplementary'))
    parser.add_argument('--extract',action='store_true',help='Also write the decompressed TSV files.')
    args=parser.parse_args();print(json.dumps(assemble_archives(args.root,args.output_dir,args.extract),indent=2))

if __name__=='__main__':main()
