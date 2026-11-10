"""Unités biologiques explicites : coordonnées natives versus alignement traduit.

BLASTX: QUERY nt, sujet aa, colonnes d'alignement aa.
TBLASTN: QUERY aa, sujet nt, colonnes aa.
TBLASTX: QUERY nt, sujet nt, colonnes aa.
NE JAMAIS assimiler 1 position nt à 1 colonne aa.
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
    def subject_label(self):return f'RÉFÉRENCE ({self.subject})'
    @property
    def alignment_label(self):return f'Colonnes alignées ({self.aligned})'


def units_for(mode:str,query_kind='ADN/ARN')->Units:
    m=mode.lower()
    if m in ('blastx','diamond-blastx'):
        return Units('nt','aa','aa','aa',True,False)
    if m=='tblastn':return Units('aa','nt','aa','aa',False,True)
    if m=='tblastx':return Units('nt','nt','aa','aa',True,True)
    if m in ('blastn','global-dna'):
        return Units('nt','nt','nt','nt',False,False)
    if m in ('blastp','diamond-blastp','mmseqs-blastp','global-aa'):
        return Units('aa','aa','aa','aa',False,False)
    unit='nt' if query_kind=='ADN/ARN' else 'aa'
    return Units(unit,unit,unit,unit,False,False)


def coord_audit(hsps,mode):
    """Audit de cohérence des coordonnées sans modification des positions BLAST.

    Calcul borné aux HSP possédant qseq et sseq ; les écarts sont diagnostiques
    (des gaps peuvent affecter le nombre de colonnes d'alignement).
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
