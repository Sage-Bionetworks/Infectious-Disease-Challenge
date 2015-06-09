from Bio import Entrez
from Bio import Geo
import pandas as pd
import synapseclient
from synapseclient import Schema, Column, Table, Row, RowSet, as_table_columns


QUERY = """("Acute respiratory distress syndrome"[All Fields] OR "communicable diseases"[MeSH Terms] OR ("communicable diseases"[MeSH Terms] OR infectious disease[All Fields]) OR ("sepsis"[MeSH Terms] OR sepsis[All Fields]) OR "sepsis"[MeSH Terms] OR "septic shock"[All Fields] OR ("pneumonia"[MeSH Terms] OR pneumonia[All Fields]) OR ("bacteremia"[MeSH Terms] OR Bacteremia[All Fields])) AND ("Homo sapiens"[Organism] AND ("50"[n_samples] : "10000000"[n_samples]))"""


syn=synapseclient.Synapse()
syn.login() 


if __name__ == '__main__':
    Entrez.email = "larsson.omberg@sagebase.org"
    handle = Entrez.esearch(db="gds", term=QUERY, retmax=500)
    record = Entrez.read(handle)
    datasets = []
    for geo_id in record['IdList']:
        handle = Entrez.esummary(db="gds", id=geo_id)
        dataset = Entrez.read(handle)[0]
        del dataset['Samples']
        del dataset['SSInfo']
        for k,v in dataset.iteritems():
            try:
                if len(v)>1000:
                    dataset[k]=v[:999]
            except TypeError:
                pass
        datasets.append(dataset)
        
df = pd.DataFrame(datasets)
df.drop(['ExtRelations', 'Projects', 'Relations'], axis=1, inplace=True)

ftpLink = Column(columnType='LINK',  maximumSize=84,  name= 'FTPLink')
columns = [ftpLink if col['name']=='FTPLink' else col for col in as_table_columns(df)]


schema = Schema(name='GEO Datasets', columns=columns, parent='syn4012977')

df.to_csv('skit.csv', encoding='utf-8', index=False)

table = syn.store(Table(schema, 'skit.csv'))



