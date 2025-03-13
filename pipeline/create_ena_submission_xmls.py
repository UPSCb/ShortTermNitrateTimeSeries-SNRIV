import xml.etree.ElementTree as ET
import pandas as pd
import hashlib
import json
import argparse
from docx import Document
import os

### PARSE JSON INPUT ### 

parser = argparse.ArgumentParser(description='Script to create XMLs file for ENA submissions')

parser.add_argument('json_file', type=str,
                    help='Input JSON file')

args = parser.parse_args()

f = open(args.json_file)

# Reads JSON file into dictionary
arg_dict = json.load(f)


### SPECIES DICTIONARY ### 

# Hard-coded dictionary with species information: taxon_id, common_name, scientific_name
species_dict = {
	"algae" : {
		"taxon_id": "2909984",
		"common_name": "Nordic green algae",
		"scientific_name": "Scenedesmus sp."
	},
	"arabidopsis": {
	  "taxon_id": "3702",
		"common_name": "thale cress",
		"scientific_name": "Arabidopsis thaliana"
	},
	"T89" : {
		"taxon_id": "47664",
		"common_name": "T89",
		"scientific_name": "Populus tremula x Populus tremuloides"
	}, 
	"yeast" : {
		"taxon_id": "4932",
		"common_name": "bakers yeast",
		"scientific_name": "Saccharomyces cerevisiae",
	},


} 

# Extract species from JSON - make sure you enter the species with its common name matching the species_dict!
species_name = arg_dict.get("species")

# Lookup the species in the nested dictionary
if species_name in species_dict:
    species_data = species_dict[species_name]
    taxon_id = species_data.get("taxon_id")
    common_name = species_data.get("common_name")
    scientific_name = species_data.get("scientific_name")

### PARSE ABSTRACT AND SAMPLE FILE ### 

# Abstract file 
abstract_file = arg_dict['abstract_file']
# Sample file
df= pd.read_csv(arg_dict['csv_file'], sep=',')


### HELPER FUNCTIONS ###

# Get unique sample names from the 'SampleName' column
def get_unique_names(file): 
   noDupes = []
   [noDupes.append(i) for i in file if not noDupes.count(i)]
   return noDupes

# Read abstract file and extract title + abstract
def extract_title_abstract(doc_path):
    expanded_path = os.path.expanduser(doc_path)
    if not os.path.isfile(expanded_path):
        raise FileNotFoundError(f"File not found: {expanded_path}")
        
    doc = Document(expanded_path)
    title = None
    abstract = None
    
    for i, para in enumerate(doc.paragraphs):
        if i == 0:  # Assuming the first paragraph is the title
            title = para.text.strip()
        elif para.text.strip().lower() == "abstract":  # Assuming "Abstract" is a separate paragraph
            abstract = doc.paragraphs[i+1].text.strip()  # The next paragraph is the abstract
    
    return title, abstract



# Read FastQ files and create md5 checksum (used if needed)
def compute_md5(file_name, block_size = 65536):
    hash_md5 = hashlib.md5()
    with open(file_name, "rb") as f:
        for chunk in iter(lambda: f.read(block_size), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

### FUNCTIONS TO CREATE XMLs ###

# Submission XML
def create_submission_xml(file): 
	root = ET.Element('SUBMISSION_SET', {"xmlns:xsi": "http://www.w3.org/2001/XMLSchema-instance", "xsi:noNamespaceSchemaLocation": "ftp://webin.ebi.ac.uk/meta/xsd/sra_1_5/SRA.submission.xsd"})
	sub = ET.SubElement(root, 'SUBMISSION', {"alias": arg_dict['alias'], 'center_name': arg_dict['center']})
	actions = ET.SubElement(sub, 'ACTIONS')
	sub_act = ET.SubElement(actions, 'ACTION')
	add_study = ET.SubElement(sub_act, 'ADD', {'source': arg_dict['alias']+".Study.xml", 'schema': 'study'})
	sub_act = ET.SubElement(actions, 'ACTION')
	add_sample = ET.SubElement(sub_act, 'ADD', {'source': arg_dict['alias']+".Sample.xml", 'schema': 'sample'})
	sub_act = ET.SubElement(actions, 'ACTION')
	add_exp = ET.SubElement(sub_act, 'ADD', {'source': arg_dict['alias']+".Experiment.xml", 'schema': 'experiment'})
	sub_act = ET.SubElement(actions, 'ACTION')
	add_run = ET.SubElement(sub_act, 'ADD', {'source': arg_dict['alias']+".Run.xml", 'schema': 'run'})
	sub_act = ET.SubElement(actions, 'ACTION')
	hold = ET.SubElement(sub_act, 'HOLD', {'HoldUntilDate': arg_dict["hold_date"]})

	tree = ET.ElementTree(root)
	ET.indent(tree, space="\t", level=0)
	tree.write(arg_dict['alias']+".Submission.xml", xml_declaration=True, encoding="utf-8")

# Sample XML
def create_sample_xml(file) : 
	num=0
	root = ET.Element('SAMPLE_SET', {"xmlns:xsi": "http://www.w3.org/2001/XMLSchema-instance", "xsi:noNamespaceSchemaLocation": "ftp://webin.ebi.ac.uk/meta/xsd/sra_1_5/SRA.sample.xsd"})
	for entry in get_unique_names(df['SampleName']):
		
		num += 1 # Increment for every sample
		length=4 # Length of sample name arg_dict['alias'] 
		sample = ET.SubElement(root, 'SAMPLE', {"alias": arg_dict['alias']+"-S" + '0' * (length - len(str(num))) + str(num), 'center_name': arg_dict['center']})
		title = ET.SubElement(sample, 'TITLE')
		title.text = entry
		sample_name = ET.SubElement(sample, 'SAMPLE_NAME')
		tax_id = ET.SubElement(sample_name, 'TAXON_ID')
		tax_id.text=taxon_id
		sci_name = ET.SubElement(sample_name, 'SCIENTIFIC_NAME')
		sci_name.text=scientific_name
		com_name = ET.SubElement(sample_name, 'COMMON_NAME')
		com_name.text=common_name
		
		desc = ET.SubElement(sample, 'DESCRIPTION')
		desc.text=df.loc[df['SampleName'] == entry, 'SampleDescription'].iloc[0]
		
		sample_attribute = ET.SubElement(sample, 'SAMPLE_ATTRIBUTES')
		attr = ET.SubElement(sample_attribute, 'SAMPLE_ATTRIBUTE')
		date=ET.SubElement(attr, 'TAG')
		date.text='collection date'
		val=ET.SubElement(attr, 'VALUE')
		val.text=arg_dict['collection_date']
		attr = ET.SubElement(sample_attribute, 'SAMPLE_ATTRIBUTE')
		loc=ET.SubElement(attr, 'TAG')
		loc.text="geographic location (country and/or sea)"
		val=ET.SubElement(attr, 'VALUE')
		val.text=arg_dict['geo_location']

	#write to file
	tree = ET.ElementTree(root)
	ET.indent(tree, space="\t", level=0)
	tree.write(arg_dict['alias']+".Sample.xml", xml_declaration=True, encoding="utf-8")

# Study XML
def create_study_xml(title, abstract, file): 
	root = ET.Element('STUDY_SET', {"xmlns:xsi":"http://www.w3.org/2001/XMLSchema-instance", "xsi:noNamespaceSchemaLocation":"ftp://webin.ebi.ac.uk/meta/xsd/sra_1_5/SRA.study.xsd"})
	study = ET.SubElement(root,'STUDY', {"alias": arg_dict['alias'], 'center_name': arg_dict['center']})
	descriptor = ET.SubElement(study, 'DESCRIPTOR')
	stitle = ET.SubElement(descriptor, 'STUDY_TITLE')
	stitle.text= title
	stype = ET.SubElement(descriptor, 'STUDY_TYPE', {'existing_study_type': arg_dict['exp_type']})
	sab = ET.SubElement(descriptor, 'STUDY_ABSTRACT')
	sab.text = abstract

	tree = ET.ElementTree(root)
	ET.indent(tree, space="\t", level=0)
	tree.write(arg_dict['alias']+".Study.xml", xml_declaration=True, encoding="utf-8")

# Experiment XML 
def create_experiment_xml(file):
	root = ET.Element('EXPERIMENT_SET', {"xmlns:xsi":"http://www.w3.org/2001/XMLSchema-instance", "xsi:noNamespaceSchemaLocation":"ftp://webin.ebi.ac.uk/meta/xsd/sra_1_5/SRA.experiment.xsd"})
	num=0
	for entry in get_unique_names(df['SampleName']):
		num += 1 # Increment for every sample
		length=4 # Length of sample name arg_dict['alias'] 
		exp = ET.SubElement(root, 'EXPERIMENT', {"alias": arg_dict['alias'] + "-E" + '0' * (length - len(str(num))) + str(num), 'center_name': arg_dict['center']})
		title = ET.SubElement(exp, 'TITLE')
		title.text = df.loc[df['SampleName'] == entry, 'ExperimentTitle'].iloc[0]
		study_ref = ET.SubElement(exp, 'STUDY_REF', {'refname': arg_dict['alias']})
		des = ET.SubElement(exp, 'DESIGN')
		des_desc = ET.SubElement(des, 'DESIGN_DESCRIPTION')
		des_desc.text = arg_dict['design_description']
		sample = ET.SubElement(des, 'SAMPLE_DESCRIPTOR', {'refname': arg_dict['alias'] + "-S" + '0' * (length - len(str(num))) + str(num)})
		lib_desc = ET.SubElement(des, 'LIBRARY_DESCRIPTOR')
		lib_name = ET.SubElement(lib_desc, 'LIBRARY_NAME')
		lib_name.text = entry
		lib_strategy = ET.SubElement(lib_desc, 'LIBRARY_STRATEGY')
		lib_strategy.text = arg_dict['lib_strategy']
		lib_source = ET.SubElement(lib_desc, 'LIBRARY_SOURCE')
		lib_source.text = arg_dict['lib_source']
		lib_select = ET.SubElement(lib_desc, 'LIBRARY_SELECTION')
		lib_select.text = arg_dict['lib_selection']
		if arg_dict["seq_read_type"] == 'PAIRED':
			lib_layout = ET.SubElement(lib_desc, 'LIBRARY_LAYOUT')
			layout = ET.SubElement(lib_layout, arg_dict["seq_read_type"], {'NOMINAL_LENGTH': arg_dict['nom_length']}) #??????
			spot_desc = ET.SubElement(des, 'SPOT_DESCRIPTOR')
			spot_desc_spec = ET.SubElement(spot_desc, 'SPOT_DECODE_SPEC')
			spot_length = ET.SubElement(spot_desc_spec, 'SPOT_LENGTH')
			spot_length.text = arg_dict['read_length']
			read_spec = ET.SubElement(spot_desc_spec, 'READ_SPEC')
			read_index = ET.SubElement(read_spec, 'READ_INDEX')
			read_index.text = '0'
			read_label = ET.SubElement(read_spec, 'READ_LABEL')
			read_label.text = 'F'
			read_class = ET.SubElement(read_spec, 'READ_CLASS')
			read_class.text = arg_dict['read_class']
			read_type = ET.SubElement(read_spec, 'READ_TYPE')
			read_type.text = 'Forward'
			base_coord = ET.SubElement(read_spec, 'BASE_COORD')
			base_coord.text = '1'
			read_spec = ET.SubElement(spot_desc_spec, 'READ_SPEC')
			read_index = ET.SubElement(read_spec, 'READ_INDEX')
			read_index.text = '1'
			read_label = ET.SubElement(read_spec, 'READ_LABEL')
			read_label.text = 'R'
			read_class = ET.SubElement(read_spec, 'READ_CLASS')
			read_class.text = arg_dict['read_class']
			read_type = ET.SubElement(read_spec, 'READ_TYPE')
			read_type.text = 'Reverse'
			base_coord = ET.SubElement(read_spec, 'BASE_COORD')
			base_coord.text = str(int(arg_dict['read_length'])+1) 
		elif arg_dict["seq_read_type"] == 'SINGLE': 
			lib_layout = ET.SubElement(lib_desc, 'LIBRARY_LAYOUT')
			spot_desc = ET.SubElement(sample, 'SPOT_DESCRIPTOR')
			spot_desc_spec = ET.SubElement(spot_desc, 'SPOT_DECODE_SPEC')
			spot_length = ET.SubElement(spot_desc_spec, 'SPOT_LENGTH')
			spot_length.text = arg_dict['read_length'] 
		platform = ET.SubElement(exp, 'PLATFORM')
		plat_type = ET.SubElement(platform, arg_dict['platform_type'])
		instrument = ET.SubElement(plat_type, 'INSTRUMENT_MODEL')
		instrument.text = arg_dict['instrument_model']
	

	#write to file
	tree = ET.ElementTree(root)
	ET.indent(tree, space="\t", level=0)
	tree.write(arg_dict['alias']+".Experiment.xml", xml_declaration=True, encoding="utf-8")

# Run XML
def create_run_xml(file): 
	root = ET.Element('RUN_SET', {"xmlns:xsi":"http://www.w3.org/2001/XMLSchema-instance", "xsi:noNamespaceSchemaLocation":"ftp://webin.ebi.ac.uk/meta/xsd/sra_1_5/SRA.run.xsd"})
	num=0
	for entry in get_unique_names(df['SampleName']):
		num += 1 # Increment for every sample
		length = 4 # Length of sample name arg_dict['alias'] 
		run = ET.SubElement(root, 'RUN', {"alias": arg_dict['alias'] + "-R" + '0' * (length - len(str(num))) + str(num), "center_name": arg_dict['center'], "run_center": arg_dict['seq_center'], "run_date": df.loc[df['SampleName'] == entry, 'SequencingDate'].iloc[0]})
		exp_ref = ET.SubElement(run, 'EXPERIMENT_REF', {"refname": arg_dict['alias'] + "-E" + '0' * (length - len(str(num))) + str(num) })
		data_block = ET.SubElement(run, 'DATA_BLOCK')
		files = ET.SubElement(data_block, 'FILES')
		# TODO: add test for md5 column
		if arg_dict["seq_read_type"] == 'PAIRED':
			if 'MD5checksum' in df.columns:
				file1 = ET.SubElement(files,'FILE',{"filename": arg_dict['alias']+"/"+df.loc[df['SampleName'] == entry, 'FileName'].iloc[0], "filetype": arg_dict['file_type'], "checksum_method": "MD5", "checksum": df.loc[df['SampleName'] == entry, 'MD5checksum'].iloc[0]})
				file2 = ET.SubElement(files,'FILE',{"filename": arg_dict['alias']+"/"+df.loc[df['SampleName'] == entry, 'FileName'].iloc[1], "filetype": arg_dict['file_type'], "checksum_method": "MD5", "checksum": df.loc[df['SampleName'] == entry, 'MD5checksum'].iloc[1]})
			else: 
				file1 = ET.SubElement(files,'FILE',{"filename": arg_dict['alias']+"/"+df.loc[df['SampleName'] == entry, 'FileName'].iloc[0], "filetype": arg_dict['file_type'], "checksum_method": "MD5", "checksum": compute_md5(df.loc[df['SampleName'] == entry, 'FileLocation'].iloc[0]+"/"+df.loc[df['SampleName'] == entry, 'FileName'].iloc[0])})
				file2 = ET.SubElement(files,'FILE',{"filename": arg_dict['alias']+"/"+df.loc[df['SampleName'] == entry, 'FileName'].iloc[1], "filetype": arg_dict['file_type'], "checksum_method": "MD5", "checksum": compute_md5(df.loc[df['SampleName'] == entry, 'FileLocation'].iloc[1]+"/"+df.loc[df['SampleName'] == entry, 'FileName'].iloc[1])})
		if arg_dict["seq_read_type"] == 'SINGLE':
			if 'MD5checksum' in df.columns:
				file = ET.SubElement(files,'FILE',{"filename": arg_dict['alias']+"/"+df.loc[df['SampleName'] == entry, 'FileName'].iloc[0], "filetype": arg_dict['file_type'], "checksum_method": "MD5", "checksum": df.loc[df['SampleName'] == entry, 'MD5checksum'].iloc[0]})
			else:
				file = ET.SubElement(files,'FILE',{"filename": arg_dict['alias']+"/"+df.loc[df['SampleName'] == entry, 'FileName'].iloc[0], "filetype": arg_dict['file_type'], "checksum_method": "MD5", "checksum": compute_md5(df.loc[df['SampleName'] == entry, 'FileLocation'].iloc[0]+"/"+df.loc[df['SampleName'] == entry, 'FileName'].iloc[0])})
		#write to file
		tree = ET.ElementTree(root)
		ET.indent(tree, space="\t", level=0)
		tree.write(arg_dict['alias']+".Run.xml", xml_declaration=True, encoding="utf-8")


# Generate XMLs based on input files 
create_submission_xml(df)
create_sample_xml(df)


title, abstract = extract_title_abstract(abstract_file)
create_study_xml(title, abstract, df)

create_experiment_xml(df)
create_run_xml(df)

