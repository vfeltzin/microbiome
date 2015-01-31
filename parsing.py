import re

def parse_otu_map_to_reads (otu_map_file):
    #Takes an otu map file and returns a dictionary with read ID as key and otu name as value.
    otu_list = list(otu_map_file)
    otu_dict = {}
    for row in otu_list:
        row = row.split()
        for item in row[1:]:
             otu_dict[item]=row[0]
    return otu_dict

def parse_otu_map (otu_map_file):
    #Takes an otu map file and returns a dictionary with otu ID as key and list of otu readIDs as value.
    otu_list = list(otu_map_file)
    otu_dict = {}
    for row in otu_list:
        row = row.split()
        otu_dict[row[0]]=row[1:]
    return otu_dict
		
def parse_tax_map (tax_map_file):
	#Takes a QIIME taxonomy map file and returns a dictionary with otu name as key and a taxonomy dictionary (which you can pass to the qiime_seq object
	tax_list = list(tax_map_file)
	tax_dict = {}
	for row in tax_list:
		tax_dict[row.split('\t')[0]] = dict(zip(['kingdom','phylum','class','order','family','genus'],row.split('\t')[1].split(';')))
	return tax_dict
	
def assign_otu_IDs_to_reads (otu_dict,qiime_seq_list):   
    for qiime_seq_obj in qiime_seq_list:
        if qiime_seq_obj.seqID in otu_dict:
            qiime_seq_obj.set_otu_ID(otu_dict[qiime_seq_obj.seqID])
    return qiime_seq_list
	
def assign_tax_IDs_to_reads (tax_dict,qiime_seq_list):
	#Takes a taxonomy map dictionary and a list of qiime_seq objects and assigns taxonomy to each qiime_seq object in the list if its otu_ID is a value in the taxonomy map dictionary
	for qiime_seq_obj in qiime_seq_list:
			if qiime_seq_obj.otu_ID in tax_dict:
				qiime_seq_obj.set_tax(tax_dict[qiime_seq_obj.otu_ID])
	return qiime_seq_list
	
def parse_qiime_fasta (qiime_fasta_file):
	#Takes a QIIME fasta file and returns a list of qiime_seq objects
	fasta_list = list(qiime_fasta_file)
	qiime_seq_list = []
	x=0
	while x<len(fasta_list):
		header = fasta_list[x].split()
		seq = fasta_list[x+1].strip('\n') 
		obj=qiime_seq(header[0].strip('>'),seq)
		obj.set_amplicon_ID(header[1])
		qiime_seq_list.append(obj)
		x = x + 2
	return qiime_seq_list
		
class qiime_seq:
	def __init__(self,seqID,seq):
		self.seqID = seqID
		self.amplicon_ID = str()
		self.otu_ID = str()
		self.seq = seq
		self.tax = {'kingdom':' ','phylum':' ','class':' ','order':' ','family':' ','genus':' ','species':' '}
		
	def set_amplicon_ID(self,amplicon_ID):
		self.amplicon_ID = amplicon_ID
		
	def set_otu_ID(self,otu_ID):
		self.otu_ID = otu_ID
		
	def set_tax(self,tax):
		self.tax['kingdom'] = tax['kingdom']
		self.tax['phylum'] = tax['phylum']
		self.tax['class'] = tax['class']
		self.tax['order'] = tax['order']
		self.tax['family'] = tax['family']
		self.tax['genus'] = tax['genus']
		#self.tax['species'] = tax['species']
	
	#To do: there should be different ways to set taxonomy on an object - you can pass a dictionary directly, assign it field by field, or take a tuple. If I think of more ways that might be useful, I'll add them later. 
	
def filter_otu_map_by_nreads(otu_map, maxreads):
	filtered_list = []
	for k,v in otu_map.items():
		if len(v)<maxreads:
			filtered_list.append(k)
			
	return(filtered_list)
	
def get_qiime_seq_obj_by_otus(otu_name_list,qiime_seq_list):
	new_qiime_seq_list = []
	for obj in qiime_seq_list:
		if obj.otu_ID in otu_name_list:
			new_qiime_seq_list.append(obj)
	return(new_qiime_seq_list)
	
def get_list_of_rep_set_IDs(rep_set_file):
	rep_set_list = list(rep_set_file)
	otu_list = []
	read_ID_list = []
	x = 0
	while x<len(rep_set_list):
		header = rep_set_list[x].split()
		otu_ID = header[0].strip('><')
		read_ID = header[1].strip('\n')
		otu_list.append(otu_ID)
		read_ID_list.append(read_ID)
		x = x + 2
	return(otu_list,read_ID_list)
	
def get_list_of_otus_from_otu_table(otu_table_file, otu_ID_prefix):
	otu_table_string = otu_table_file.read()
	regexp = otu_ID_prefix+'\d+'
	otu_ID_list = re.findall(regexp,otu_table_string)
	return(otu_ID_list)
	
def parse_qiime_fasta_selectively(read_IDs,qiime_fasta_file):
	fasta_list = list(qiime_fasta_file)
	qiime_seq_list = []
	x=0
	while x<len(fasta_list):
		header = fasta_list[x].split()
		seq = fasta_list[x+1].strip('\n')
		if header[0].strip('>') in read_IDs:
			obj=qiime_seq(header[0].strip('>'),seq)
			obj.set_amplicon_ID(header[1])
			qiime_seq_list.append(obj)
		x = x + 2
	return(qiime_seq_list)
	
def parse_qiime_fasta_from_repset(qiime_fasta_file):
    fasta_list = list(qiime_fasta_file)
    qiime_seq_list = []
    x=0
    while x<len(fasta_list):
        header = fasta_list[x].split()
        seq = fasta_list[x+1].strip('\n')
        obj=qiime_seq(header[1],seq)
        obj.set_otu_ID(header[0].strip('><'))
        qiime_seq_list.append(obj)
        x=x+2
    return(qiime_seq_list)
		
		
