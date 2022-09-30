import untangle
import sys
import pandas as pd
import re

xml_file = sys.argv[1]
#xml_file = 'combined.xml'
# obj = untangle.parse(xml_file)
#
# biosample_attributes = list(obj.DocumentSummarySet.DocumentSummary.SampleData.BioSample.Attributes.Attribute)
# attribute_dict = {a.get_attribute('attribute_name'): a.cdata for a in biosample_attributes}
# organism = obj.DocumentSummarySet.DocumentSummary.Organism.cdata
# attribute_dict['organism'] = organism
# print(attribute_dict)

# For reading a whole xml file (new method)
with(open(xml_file, 'r') as f):
    newItemBool = True
    xml_list = []
    xml_index = 0
    for i, line in enumerate(f.readlines()):
        if i>3: # skip header...
            #print(line)
            if line.strip()=="</DocumentSummarySet>":
                print('end')
                xml_index = 'end'
            if line.strip()=="<DocumentSummary>":
                xml_list.append('')
            if xml_index!='end':
                xml_list[xml_index] += line
            if line.strip()=="</DocumentSummary>":
                xml_index += 1


biosample_dict = {}
for i, xml_string in enumerate(xml_list):
    obj = untangle.parse(xml_string)
    biosample_accession = re.sub('BioSample:', '', str(obj.DocumentSummary.SourceSample.cdata))
    biosample_attributes = list(obj.DocumentSummary.SampleData.BioSample.Attributes.Attribute)
    attribute_dict = {a.get_attribute('attribute_name'): a.cdata for a in biosample_attributes}
    organism = obj.DocumentSummary.Organism.cdata
    attribute_dict['organism'] = organism
    biosample_dict[biosample_accession] = attribute_dict
    #print(i, biosample_accession, attribute_dict)

#
useful_columns = ['organism', 'host', 'collected_by',  'geo_loc_name', 'collection_date', 'lat_lon']

# 'ENA first public', 'INSDC first public'

df = pd.DataFrame.from_dict(biosample_dict, orient='index')

df[useful_columns].to_csv(xml_file+'.csv')
# create object for each of the separated files - read from <DocumentSummary> to </DocumentSummary>
