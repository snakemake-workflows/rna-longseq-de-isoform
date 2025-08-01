import sys
import csv
import requests
import time

# Start logging
sys.stderr = sys.stdout = open(snakemake.log[0], "w")


input_file = snakemake.input[0]
output_file = snakemake.output[0]


# Obtain the protein name from the UniRef databank ID
def get_protein_name(uniref_id):
    url = f"https://rest.uniprot.org/uniref/{uniref_id}"
    # Check if the url is available
    try:
        response = requests.get(url, timeout=30)
        # If available get the protein name for the UniRef ID
        if response.status_code == 200:
            data = response.json()
            rep_member = data.get("representativeMember", {})
            protein_name = rep_member.get("proteinName", "Unknown")
            return protein_name
        # If the url is unavailable throw error
        else:
            print(f"Failed for {uniref_id} with status {response.status_code}")
            return "Request failed"
    except Exception as e:
        print(f"Error: {e}")
        return "Error"


# Creates a csv file containing: gene name from the annotation, UniRef ID, protein name
def process_blast_output(input, output):
    seen_ids = set()
    # Opens the lambda results input and writes a header for the output file
    with open(input, "r") as infile, open(output, "w", newline="") as outfile:
        reader = csv.reader(infile, delimiter="\t")
        writer = csv.writer(outfile)
        writer.writerow(["Gene", "UniRefID", "ProteinName"])

        # Iterates through each row of the input and uses the first two entries to get gene and uniref_id
        for row in reader:
            if len(row) < 2:
                continue
            gene, uniref_id = row[0], row[1]
            if uniref_id in seen_ids:
                continue
            seen_ids.add(uniref_id)

            # Uses get_protein-name to obtain and write the protein name
            protein_name = get_protein_name(uniref_id)
            writer.writerow([gene, uniref_id, protein_name])
            time.sleep(1)


# Starts the annotation
process_blast_output(input_file, output_file)
