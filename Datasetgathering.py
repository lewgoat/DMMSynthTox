import pubchempy as pcp
import csv
import requests


def get_data(cid):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/json/"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        data = str(data)
        dataout ={'cid':cid, 'ld50':'' , 'smiles':''} #create a dictionary to store the data
        #find the index of the ld50 in data string
        try:
            index = data.find("LD50 Rat")
            dataout['ld50'] = data[index+14:index+24]
        except:
            dataout['ld50'] = None
    else:
        #return a request error

        return TimeoutError("Error: API request unsuccessful.")


# an empty list to store the compound ids
compound_ids = []
with open('PubChem_compound_list_aS7My_0hmJ2vsxCqktJZgGjZFLlr-0oYMD1RVCssQ1UrNX8.csv', newline='') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        compound_ids.append(row[0])
#remove the first element of the compound ids
compound_ids.pop(0)
#remove the leading empty spcae of the first element
compound_ids[0] = compound_ids[0].lstrip()

datatable = []

#get and print the first 10 ld50 values
for i in range(10):
    datatable.append(get_data(compound_ids[i]))
    
print(datatable)
