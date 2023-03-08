import pubchempy as pcp
import csv
import requests
count = 0
#multi threaded function to get the ld50 value from the pubchem api


def get_data(cid):
    print(str(count)+"/108141")
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/json/"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        data = str(data)
        #find the index of the ld50 in data string
        try:
            index = data.find("LD50 Rat")
            
            dataout = data[index:index+45]
            #remove commas from the ld50 value
            dataout = dataout.replace(',', '')
            return dataout
        except:
            return None
    else:
        #return a request error

        return TimeoutError("Error: API request unsuccessful.")


#function that removes unwanted characters from the ld50 value
def remove_chars(string):
    '''
    removes unwanted characters from the ld50 value
    '''
    try:
        string = string.replace('LD50 Rat', '')
        string = string.replace('(', '')
        string = string.replace(')', '')
        string = string.replace('{', '')
        string = string.replace('}', '')
        string = string.replace('[', '')
        string = string.replace(']', '')
        string = string.replace(' ', '')
        string = string.replace(',', '')
        string = string.replace('"', '')
        string = string.replace('oral', '')
        #remove any non numeric characters that are not m,u,g or k or -
        
        string = ''.join([i for i in string if i.isdigit() or i == '.' or i == 'm' or i == 'u' or i == 'g' or i == 'k' or i == '-'])
        
        #an edge case condition where there is a range given for the ld50 value
        if string.find('-') != -1:
            #split the string at the - and keep the first element and nonnumeric characters in second element
            string = string.split('-')[0] + ''.join([i for i in string.split('-')[1] if i.isalpha()])
        
        print(string)
        #keep up to second g
        string = string[:string.find('g')+4]
        return string
    except:
        return None
def consistent_units(ldval):
    '''
    converts the units of the ld50 value to g/kg and removes the units
    '''
    try:
    #if it contains more than one . remove all but the last one
        if ldval.count('.') > 1:
            #find the index of the last .
            index = ldval.rfind('.')
            #remove all full stops before the last full stop
            ldval = ldval[:index] + ldval[index:].replace('.', '')                     


        if ldval.find('mg') != -1:
            #convert mg to g
            
            ldval = ldval.replace('mg', '')
            ldval = ''.join([i for i in ldval if i.isdigit() or i == '.'])
            ldval = float(ldval) / 1000
            
            return ldval
        elif ldval.find('ug') != -1:
            #convert ug to g
            ldval = ldval.replace('ug', '')
            #remove non-numeric characters
            ldval = ''.join([i for i in ldval if i.isdigit() or i == '.'])
            ldval = float(ldval) / 1000000
            return ldval
        elif ldval.find('g/kg') != -1:
            #remove g/kg from the string
            ldval = ''.join([i for i in ldval if i.isdigit() or i == '.'])
            ldval = ldval.replace('g/kg', '')
            ldval = float(ldval)
            return ldval
    except:
        return None
    


#open cid2smiles text file and convert to dictionary
 


# an empty list to store the compound ids
compound_ids = []
with open('FullListCIDpreprocess.csv', newline='') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        compound_ids.append(row[0])
#remove the first element of the compound ids
compound_ids.pop(0)
#remove the leading empty spcae of the first element
compound_ids[0] = compound_ids[0].lstrip()
#save compound ids to a csv file
with open('compound_ids.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(compound_ids) 
datatable = []

cid2smiles = {}
with open('cid2smiles.txt', 'r') as f:
    for line in f:
        (key, val) = line.split()
        cid2smiles[int(key)] = val
#preallocate the datatable list and fill it with the cid and smiles
datatable= [{'cid':cid, 'ld50':'' , 'smiles':cid2smiles[int(cid)]} for cid in compound_ids]
#loop through the compound ids and get the ld50 value using the get_data function
for i in range(len(compound_ids)):
    datatable[i]['ld50'] = get_data(compound_ids[i])
    count += 1

#save datatable to file
with open('datatablePREPROCESSED.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(datatable)
cid2smiles = {}
#open datatable file and convert to dictionary
datatable = []
with open('datatablePREPROCESSED.csv', 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        datatable.append(row)

#remove unwanted characters from the ld50 value using remove chars function
for i in range(len(datatable)):
    datatable[i]['ld50'] = remove_chars(datatable[i]['ld50'])
    if datatable[i]['ld50'] == 'None':
        datatable[i]['ld50'] = None

#conistent units
for i in range(len(datatable)):
    if datatable[i]['ld50'] != None:
        datatable[i]['ld50'] = consistent_units(datatable[i]['ld50'])
    else:
        datatable[i]['ld50'] = None
#remove any compounds that have no ld50 value
datatable = [x for x in datatable if x['ld50'] != None]



#convert datatable from dictionary to string with the keys as the first row
datatable = [list(datatable[0].keys())] + [list(x.values()) for x in datatable]
#convert datatable to csv format 
datatable = [','.join(map(str, x)) for x in datatable]
# add a new line to the end of each element in the datatable list
datatable = [x + '\n' for x in datatable]
#convert datatable to a string and remove double quotes
datatable = ''.join(datatable)
datatable = datatable.replace('"', '')

print(datatable)


#save datatable to a csv file
with open('datatable.csv', 'w') as f:
    f.write(datatable)