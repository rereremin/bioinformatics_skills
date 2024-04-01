import pandas as pd
from Bio import Entrez

def fetch_annotation(gene_id):
    try:
        Entrez.email = "your.email@example.com"  # Enter your email address here
        handle = Entrez.efetch(db="gene", id=gene_id, rettype="gb", retmode="text")
        record = handle.read()
        handle.close()
        return record
    except Exception as e:
        print(f"Error fetching annotation for gene ID {gene_id}: {e}")
        return None

def main():
    with open("gene_ids.txt", "r") as file:
        gene_ids = [line.strip() for line in file if line.strip()]

    data = []
    for gene_id in gene_ids:
        print(f"Fetching annotation for gene ID: {gene_id}")
        annotation = fetch_annotation(gene_id)
        if annotation:
            annotation = annotation.split('\n')[2]
            data.append({"Gene ID": gene_id, "Annotation": annotation})
    df = pd.DataFrame(data)
    return df

df = main()
df.to_csv('gene_annotation.csv')
