from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import os

def qiimeseq_blaster (qiime_seq_list):
    "Takes a list of qiime_seq objects and returns a list of blast_record objects."
    blast_record_list = []
    for obj in qiime_seq_list:
        result_handle = NCBIWWW.qblast("blastn", "nt", obj.seq)
        blast_record = NCBIXML.read(result_handle)
        blast_record_list.append(blast_record)
        
    return(blast_record_list)
		