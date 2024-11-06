import pandas as pd


def dissociations_table_insert(diss_method, conn):
    # check if the record already exists
    cursor = conn.cursor()
    sql = "SELECT * FROM dissociations WHERE dissociation_method = ?"
    val = (str(diss_method), )
    cursor.execute(sql, val)
    res = cursor.fetchone()
    if res is None:
        sql = "INSERT INTO dissociations (dissociation_method) VALUES (?)"
        val = (str(diss_method), )
        cursor.execute(sql, val)
        conn.commit()
        print("Dissociation record added!")
    else:
        print("Dissociation record already exists")
        dissociation_id = res[0]
        print(f"current dissociation id is : {dissociation_id}")
    
    
# insert a vendor record to vendors table
def vendors_table_insert(vendor_name, conn):
    # check if the record already exists
    cursor = conn.cursor()
    sql = "SELECT * FROM vendors WHERE vendor_name = ?"
    val = (str(vendor_name), )
    cursor.execute(sql, val)
    res = cursor.fetchone()
    if res is None:
        sql = "INSERT INTO vendors (vendor_name) VALUES (?)"
        val = (str(vendor_name), )
        cursor.execute(sql, val)
        conn.commit()
        print("Vendor record added!")
        sql = "SELECT vendor_id FROM vendors WHERE vendor_name = ?"
        val = (str(vendor_name), )
        cursor.execute(sql, val)
        vendor_res = cursor.fetchone()
        vendor_id = vendor_res[0]
    else:
        print("Vendor record already exists")
        vendor_id = res[0]
        print(f"current vendor id is : {vendor_id}")
    return vendor_id
    
    
# insert an instrument record to instruments table
def instruments_table_insert(vendor_id, instrument_name, conn):
    # check if the record already exists
    cursor = conn.cursor()
    sql = "SELECT * FROM instruments WHERE vendor_id= ? AND instrument_name = ?"
    val = (int(vendor_id), str(instrument_name))
    cursor.execute(sql, val)
    res = cursor.fetchone()
    if res is None:
        sql = "INSERT INTO instruments (vendor_id, instrument_name) VALUES (?, ?)"
        val = (int(vendor_id), str(instrument_name))
        cursor.execute(sql, val)
        conn.commit()
        print("Instrument record added!")
    else:
        print("Instrument record already exists")
        instrument_id = res[0]
        print(f"current instrument id is : {instrument_id}")

    
# insert a species record to species table
def species_table_insert(species_name, taxon_id, conn):
    cursor = conn.cursor()
    # check if the record already exists
    sql = "SELECT * FROM species WHERE species_name = ?"
    val = (species_name, )
    cursor.execute(sql, val)
    res = cursor.fetchone()
    if res is None:
        sql = "INSERT INTO species (species_name, taxon_id) VALUES (?, ?)"
        val = (str(species_name), int(taxon_id))
        cursor.execute(sql, val)
        conn.commit()
        print("Species record added!")
    else:
        print("Species record already exists")
        species_id = res[0]
        print(f"current species id is : {species_id}")
   
    
def species_diss_instru_vendor_tables(lib_parameters, conn):
    vendor_lst = lib_parameters['vendor_name']
    dissociation_lst = lib_parameters['dissociation_method']
    instrument_lst = lib_parameters['instrument_name']  
    species_taxonID_lst = lib_parameters['species_taxon_id']        
    # generate species table
    species_data = pd.DataFrame(species_taxonID_lst)
    for i in range(len(species_data)):
        species_name, taxon_id = species_data.loc[i, ['common_name','taxon_id']]
        species_table_insert(species_name, taxon_id, conn)
    # generate vendor, instrument table
    for i in range(len(vendor_lst)):
        vendor_name = vendor_lst[i]
        vendor_id = vendors_table_insert(vendor_name, conn)
        for j in range(len(instrument_lst)):
            instrument_name = instrument_lst[j]
            instruments_table_insert(vendor_id, instrument_name, conn)
    # generate dissociation table
    for i in range(len(dissociation_lst)):
        dissociation_method = dissociation_lst[i]
        dissociations_table_insert(dissociation_method, conn)
             

