import pandas as pd
import sys
def merge_tables(tab1, tab2, otab):    
    tab1= pd.read_csv(tab1, sep= "\t", low_memory= False, on_bad_lines= "skip")
    tab2= pd.read_csv(tab2, sep= "\t", low_memory= False, on_bad_lines= "skip")
    mtab= pd.merge(tab1, tab2, how= "left", on= ["CHROM", "POS", "ID", "REF", "ALT"])
    new_cols= { "CHROM": "CHROMOSOME", "POS": "POSITION", "ID": "ID", "REF": "REFERENCE_ALLELE", "ALT": "ALTERNATE_ALLELE",
                "ANN[*].EFFECT": "EFFECT", "ANN[*].IMPACT": "IMPACT", "ANN[*].GENE": "GENE", 
                "ANN[*].GENEID": "GENEID", "ANN[*].FEATURE": "FEATURE", 
                "ANN[*].FEATUREID": "FEATUREID", "ANN[*].HGVS_C": "CHROMOSOME_HGVS", "ANN[*].HGVS_P": "PROTEIN_HGVS",
                "CLNSIG": "CLINICAL_SIGNIFICANCE", "CLNHGVS": "CLINVAR_HGVS", "CLNDN": "CLINVAR_DISEASE_NAME"}
    mtab= mtab.fillna("-")
    mtab.rename(columns= new_cols, inplace= True)
    mtab.to_csv(otab, sep= "\t", index= False)

tab1= sys.argv[1]
tab2= sys.argv[2]
otab= sys.argv[3]
merge_tables(tab1, tab2, otab)