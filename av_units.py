"""Explicit biological units: native coordinates versus translated alignment.

BLASTX: QUERY nt, reference aa, alignment columns aa.
TBLASTN: QUERY aa, sujet nt, columns aa.
TBLASTX: QUERY nt, sujet nt, columns aa.
Never equate one nucleotide position with one amino-acid column.
"""
from __future__ import annotations
from dataclasses import dataclass

@dataclass(frozen=True)
class Units:
    query: str
    subject: str
    aligned: str
    msa: str
    translated_query: bool
    translated_subject: bool
    @property
    def query_label(self):return f'QUERY ({self.query})'
    @property
    def subject_label(self):return f'REFERENCE ({self.subject})'
    @property
    def alignment_label(self):return f'Aligned columns ({self.aligned})'


def units_for(mode:str,query_kind='DNA/RNA')->Units:
    m=mode.lower()
    if m in ('blastx','diamond-blastx'):
        return Units('nt','aa','aa','aa',True,False)
    if m=='tblastn':return Units('aa','nt','aa','aa',False,True)
    if m=='tblastx':return Units('nt','nt','aa','aa',True,True)
    if m in ('blastn','global-dna'):
        return Units('nt','nt','nt','nt',False,False)
    if m in ('blastp','diamond-blastp','mmseqs-blastp','global-aa'):
        return Units('aa','aa','aa','aa',False,False)
    unit='nt' if query_kind=='DNA/RNA' else 'aa'
    return Units(unit,unit,unit,unit,False,False)


def coord_audit(hsps,mode):
    """Check coordinate consistency without changing BLAST positions.

    Limited to HSPs with qseq and sseq; differences are diagnostic
    (gaps can affect the number of alignment columns).
    """
    units=units_for(mode)
    rows=[]
    for h in hsps:
        qlen=h.qspan[1]-h.qspan[0]+1
        slen=h.sspan[1]-h.sspan[0]+1
        qnongap=sum(c!='-' for c in h.qseq)
        snongap=sum(c!='-' for c in h.sseq)
        qfactor=3 if units.translated_query else 1
        sfactor=3 if units.translated_subject else 1
        rows.append(dict(subject=h.subject,qstart=h.qstart,qend=h.qend,sstart=h.sstart,send=h.send,
            query_unit=units.query,reference_unit=units.subject,alignment_unit=units.aligned,
            query_span=qlen,reference_span=slen,query_aligned_chars=qnongap,reference_aligned_chars=snongap,
            expected_query_span=qnongap*qfactor if h.qseq else None,
            expected_reference_span=snongap*sfactor if h.sseq else None,
            query_consistent=abs(qlen-qnongap*qfactor)<=2 if h.qseq else None,
            reference_consistent=abs(slen-snongap*sfactor)<=2 if h.sseq else None,
            qframe=h.qframe,sframe=h.sframe))
    return rows
