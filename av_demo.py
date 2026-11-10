"""Bundled examples available in wheels and portable zipapps."""
from pathlib import Path

SAMPLES = {'example_query.fna': '>query_ASV_test\nATGAAAGTTGACGACTGGCGCATGAGATTGGATGACAAGCTGGGTGCTACCGTTGCTATGGCCGCATCCCGTACTAGCAGACTGTGGCGTCGGAGACGTTAGCTTCTTCGACATGATCAAGATCCGCCAGGCCGACTTCGCGAAGCTTGTGGGAGATACCGTAGCGTGCCTGATGGACGCTTTCTTCGACGGATCAGTGGAGAAGGCTCCGATGGCAGAACCTCGTGGACTCAGCAAATCCATGCTGCCTGACGCGTCCGAATGA\n', 'example_msa.fasta': '>seq_1\nATGCA-TACCG\n>seq_2\nATGCG-TACCG\n>seq_3\nATGCGTTACGG\n>seq_4\nATGCA-TACGG\n', 'example_refs.faa': '>ref_HSP_1 Example homologous protein\nMPLDIGDAVKRHITATDGKVLWNIA\n>ref_HSP_2 second protein\nMPLDIGDAVKRHITATDGRVLWNIA\n>ref_HSP_3 divergent protein\nMPLDIADAVKRHITATDGRVLWKIA\n'}

def write_demo(folder):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    for name,content in SAMPLES.items():(folder/name).write_text(content,encoding='utf-8')
    return folder/'example_msa.fasta'
