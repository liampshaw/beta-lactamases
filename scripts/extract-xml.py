import untangle
import sys

xml_file = sys.argv[1]

obj = untangle.parse(xml_file)

biosample_attributes = list(obj.DocumentSummarySet.DocumentSummary.SampleData.BioSample.Attributes.Attribute)
attribute_dict = {a.get_attribute('attribute_name'): a.cdata for a in biosample_attributes}
organism = obj.DocumentSummarySet.DocumentSummary.Organism.cdata
attribute_dict['organism'] = organism
print(attribute_dict)

# For reading a whole xml file (new method)
with(open(xml_file, 'r') as f):
    for i, line in enumerate(f.readlines()):
        if i>3: # skip header...
            if line.strip()=="<DocumentSummary>":
                newItemBool = True
            # create object for each of the separated files - read from <DocumentSummary> to </DocumentSummary>
   
biosample_attributes = list(obj.DocumentSummary.SampleData.BioSample.Attributes.Attribute)
attribute_dict = {a.get_attribute('attribute_name'): a.cdata for a in biosample_attributes}
organism = obj.DocumentSummary.Organism.cdata
attribute_dict['organism'] = organism

       
