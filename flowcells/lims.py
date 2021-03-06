#!/bin/env python

import sys
import re
import argparse
import pandas as pd
import pygsheets
import os
import glob


#Functions for accessing Maurano Lab LIMS through pygsheets


#NB Need to get 
#https://pygsheets.readthedocs.io/en/latest/authorization.html
#create project, service account key, share sheet with service account manually from google drive UI


#https://stackoverflow.com/questions/26987222/checking-whitespace-in-a-string-python/26987329
import string
def contains_whitespace(s):
    return True in [c in s for c in string.whitespace]


#Convenience function
def getValueFromLIMS(lims, bs, colname):
    if lims is None:
        #Shouldn't happen but fail gracefully anyway
        print("Couldn't get ", colname, " from LIMS!", sep="", file=sys.stderr)
        return ""
    else:
        matches = lims[lims['Sample #'] == bs][colname].unique()
        if len(matches) == 0:
            #Multiple matches should never occur if BS numbers are unique
            raise Exception("Found no LIMS entries for " + bs)
        if len(matches) > 1:
            #Multiple matches should never occur if BS numbers are unique
            raise Exception("Found multiple LIMS entries for " + bs)
#        print('LIMS: ', bs, ":", colname, '=>', matches[0], ".", sep="")
        return matches[0]

#Pull the LIMS sheet from google using the service account secrets file and spreadsheet ID.
def getLIMSsheet(sheet):
    try:
        #Read google sheets ID from the file system
        lims_gsheets_ID_file = open('/vol/mauranolab/flowcells/LIMS.gsheets_ID.txt', 'r')
        lims_gsheets_ID = lims_gsheets_ID_file.readlines()[0].rstrip("\n")
        lims_gsheets_ID_file.close()
        
        gc = pygsheets.authorize(service_file='/vol/mauranolab/flowcells/mauranolab.json')
        sh = gc.open_by_key(lims_gsheets_ID)
        wks = sh.worksheet_by_title(sheet)
        df = wks.get_as_df(value_render="UNFORMATTED_VALUE")
        df.fillna(value='', inplace=True) #Maybe pandas v1 problem? Or just change '' to None?
        #Mask per-FC headers and space between FCs (any row that is only empty lines)
        mask = ~df['Sample Name'].str.startswith('#') & (df!="").any(1)
    except Exception as e:
        #Doesn't print exception name right, e.g. if header is corrupted in google docs: "AttributeError: 'KeyError' object has no attribute 'argument'"
        print("WARNING could not load sheet " + sheet + " from google sheets: ", e.message, '\n', e.argument)
        wks = None
        df = None
    return wks, df, mask


#Verify consistency in common entries between Sample Sheet and LIMS sheets
#Projects argument suppresses less important inconsistencies unless project is on that comma-separated list
#TODO enforce illegal characters in sample name?
def validateSampleSheetAgainstLIMS(lims, seq, limsMask, seqMask, projects="Maurano,CEGS"):
    projectList = projects.split(",")
    
    print("Level", "Description", "Sample Name", "Sample #", "Key", "LIMS_value", "SampleSheet_value", sep="\t")
    
    print("#Check sequencing sheet for consistency with LIMS. Projects=", projects, sep="")
    commonCols = set(lims.columns.values).intersection(set(seq.columns.values))
    numMissingSamples = 0
    numMultipleSamples = 0
    #Iterate through the sequencing sheet one row at a time, pulling any matching rows from LIMS
    #This validation only applies to samples that also appear on the sequencing sheet
    for seqRow in seq.index.values[seqMask]:
        curSeq = seq.iloc[seqRow]
        bs = curSeq['Sample #']
        SampleName = curSeq['Sample Name']
        
        curLims = lims[limsMask & lims['Sample #'].isin([bs])]
        numEntriesInLIMS = curLims.shape[0]
        if numEntriesInLIMS == 0:
            print("ERROR", "Can't find row " + str(seqRow) + " in LIMS!", SampleName, bs, "", "", "", sep="\t")
            numMissingSamples += 1
        elif numEntriesInLIMS >= 2:
            print("ERROR", "found " + str(numEntriesInLIMS) + " entries in LIMS!", SampleName, bs, "", "", "", sep="\t")
            numMultipleSamples += 1
        else:
            #Exactly one LIMS entry, but curLims still needs to be accessed using .values.item():
            #Only check additional metadata for specified projects to avoid excess verbiage
            if projects=='' or curLims['Lab'].values.item() in projectList:
                for col in commonCols:
                    if curSeq[col] != curLims[col].values.item():
                        if curSeq[col] == "" or curLims[col].values.item() == "":
                            print("WARNING", "missing info", SampleName, bs, col, curLims[col].values.item(), curSeq[col], sep="\t")
                        else:
                            print("ERROR", "inconsistent info", SampleName, bs, col, curLims[col].values.item(), curSeq[col], sep="\t")
                
                for col in set(seq.columns.values):
                    if isinstance(curSeq[col], str) and curSeq[col] != curSeq[col].strip():
                        print("WARNING", "leading/trailing whitespace in sequencing sheet", SampleName, bs, col, "", curSeq[col], sep="\t")
    
    numDuplicateSamples = 0
    for index, curLims in lims[lims.duplicated(subset="Sample #")].iterrows():
        bs = curLims['Sample #']
        SampleName = curLims['Sample Name']
        print("ERROR", "duplicate BS number", SampleName, bs, "", "", "", sep="\t")
        numDuplicateSamples += 1
    
    print("#Check LIMS for integrity. Projects=", projects, sep="")
    #Iterate through LIMS
    #This validation applies to all samples in LIMS
    lastBS = None
    for limsRow in lims.index.values[limsMask]:
        curLims = lims.iloc[limsRow]
        bs = curLims['Sample #']
        SampleName = curLims['Sample Name']
        
        if lastBS is not None and bs <= lastBS:
            print("ERROR", "sample out of order", SampleName, bs, "Sample #", "", "", sep="\t")
        lastBS=bs
        
        #Only check additional metadata for specified projects to avoid excess verbiage
        if projects=='' or curLims['Lab'] in projectList:
            for col in set(lims.columns.values):
                if str(curLims[col]) != str(curLims[col]).strip():
                    print("WARNING", "leading/trailing whitespace in LIMS", SampleName, bs, col, curLims[col], "", sep="\t")
            
            for col in ['Genetic Modification', 'Bait set', 'Custom Reference']:
                if contains_whitespace(str(curLims[col])):
                    print("WARNING", "internal whitespace", SampleName, bs, col, curLims[col], "", sep="\t")
            
            for col in ["Parent Library", "Pool ID"]:
                if curLims[col] != "":
                    #Require exactly 1 match in LIMS
                    if lims[limsMask & lims['Sample #'].isin([curLims[col]])].shape[0] != 1:
                        print("ERROR", "invalid " + col, SampleName, bs, col, curLims[col], "", sep="\t")
            
            if curLims["Genetic Modification"] != "" and curLims["Custom Reference"] != "":
                geneticModifications = set([ re.sub(r"\[.+\]$", "", cur) for cur in curLims["Genetic Modification"].split(",") ])
                customReferences = set(curLims["Custom Reference"].split(","))
                if geneticModifications & customReferences:
                    print("ERROR", "duplicate custom reference", SampleName, bs, "Genetic Modification/Custom Reference", geneticModifications & customReferences, "", sep="\t")
            
            requiredColsBySampleType = { "DNA Capture": ["Parent Library", "Bait set", "Pool ID"] }
            for sampleType in requiredColsBySampleType:
                for col in requiredColsBySampleType[sampleType]:
                    if curLims["Sample Type"] in [sampleType]:
                        if curLims[col] == "":
                            print("WARNING", col + " column is required for sample type " + sampleType, SampleName, bs, col, curLims[col], "", sep="\t")
            
            colsSpecificToSampleType = { "R1 Trim (P5)": ["Amplicon", "Transposon DNA", "Transposon RNA", "Transposon iPCR", "Transposon 10xRNA"], "R2 Trim (P7)": ["Amplicon", "Transposon DNA", "Transposon RNA", "Transposon iPCR", "Transposon 10xRNA"], "Bait set" : ["DNA Capture", "Amplicon"] }
            for col in colsSpecificToSampleType:
                if curLims[col] != "":
                    if curLims["Sample Type"] not in colsSpecificToSampleType[col]:
                        print("WARNING", col + " column is not allowed for sample type " + curLims["Sample Type"], SampleName, bs, col, curLims[col], "", sep="\t")
    
    
    print("", file=sys.stderr)
    print("#Clone info")
    #Leave out "Species" right now since it is mainly for mappings
    clones = lims[["Clone ID", "Genetic background / Individual", "Sex", "Genetic Modification"]].drop_duplicates()
    clones = clones[clones["Clone ID"] != ""]
    for index, curClones in clones[clones.duplicated(subset="Clone ID")].iterrows():
        cloneid = curClones['Clone ID']
        print("", file=sys.stderr)
        print("ERROR", "Clone with inconsistent information", cloneid, sep="\t", file=sys.stderr)
        print(clones[clones["Clone ID"] == cloneid], file=sys.stderr)
    
    
    #Print summary at end so it isn't missed
    print("", sep="", file=sys.stderr)
    print(str(numMissingSamples) + " total samples missing from LIMS sheet", sep="", file=sys.stderr)
    print(str(numMultipleSamples) + " total samples matching multiple entries in LIMS sheet", sep="", file=sys.stderr)
    print(str(numDuplicateSamples) + " total samples with duplicate BS numbers", sep="", file=sys.stderr)


#Pass in a data frame containing updates to make 
#Rows will be matched by 'Sample #', data in remaining columns will be transfered to the live copy of LIMS or Sample Sheet
#df must be the panda df mirror of wks, used to get coords efficiently
#TODO safer to make changes offline then sync? Not clear if this is implemented yet, unlink() seems to break get_value
#I think you do:
#wks.unlink()
#...
#wks.link(syncToCloud=True)
def updateSheetFromTable(wks, df, updates, commit=True):
    if any(updates['Sample #'].duplicated()):
        print("ERROR: found duplicate Sample #s in update table!")
        return
    errors = 0
    colnames = df.columns.values
    excludedCols = ['Sample #'] #Don't update data in these columns
    for BS in updates['Sample #']:
        if re.match("^BS[0-9]+[A-Z]$", BS) is None:
            print("WARNING: empty BS provided" + BS)
            continue
        #kind of a hassle to get column access
        #loc = wks.find(BS, matchEntireCell=True)
        rows = df.index[df['Sample #'] == BS].values
        #I think row + 2 because of 0 -> 1 based indexing plus the header row
        for row in rows + 2:
            print("\nRow " + str(row) + ":" + BS)
            curBS = str(wks.get_value((row, df.columns.get_loc('Sample #')+1)))
            if BS != curBS:
                errors += 1
                print("ERROR: " + curBS + " does not match expected sample " + BS)
            else:
                for col in set(colnames).intersection(set(updates.columns.values)).difference(excludedCols):
                    coords = (row, df.columns.get_loc(col)+1)
                    oldvalue = wks.get_value(coords)
                    newvalue = updates[updates['Sample #'] == BS][col].values.item()
                    if oldvalue != newvalue:
                        print(col + "=" + str(oldvalue) + " => " + str(newvalue))
                        if commit:
                            wks.update_value(coords, newvalue)
    if errors > 0:
        print("ERROR: "+ errors + " found!")


#Search and replace functionality for FC info
#BUGBUG I don't think this is safe to run twice in a given session
def replaceFCinfo(seqWks, seq, key, value, newvalue, commit=False):
    fcMask = seq['Sample Name'].str.startswith('#')
    curFC = None
    for seqRow in seq.index.values[fcMask]:
        curSeq = seq.iloc[seqRow]
        if curSeq['Sample Name']=="#Barcode":
            #BUGBUG FC is customarily after run name, so it won't be set for the first line
            curFC = curSeq['Sample #']
        elif curSeq['Sample Name']=="#"+key:
            #I think row + 2 because of 0 -> 1 based indexing plus the header row
            coords = (seqRow+2, seq.columns.get_loc('Sample #')+1)
            oldvalue = seqWks.get_value(coords)
            if (value is None or oldvalue == value) and oldvalue!=newvalue:
                print("Row ", seqRow+2, " (", curFC, "): ", key, "=", str(oldvalue), " => ", str(newvalue), sep="")
                if commit:
                    seqWks.update_value(coords, newvalue)


#Remove empty FC info records matching a given key
def removeEmptyFCinfo(seqWks, seq, key, commit=True):
    fcMask = seq['Sample Name'].str.startswith('#')
    #Iterate in reverse order since updating row index to account for deleted rows does not seem to work
    for seqRow in reversed(seq.index.values[fcMask]):
        curSeq = seq.iloc[seqRow]
        if curSeq['Sample Name']=="#"+key and curSeq['Sample #']=="":
            print("Removing row ", seqRow, ": ", curSeq['Sample Name'], sep="")
            if commit:
                #I think row + 2 because of 0 -> 1 based indexing plus the header row
                seqWks.delete_rows(int(seqRow+2), 1)


###Command line operation
if __name__ == "__main__":
    limsWks, lims, limsMask = getLIMSsheet("LIMS")
    seqWks, seq, seqMask = getLIMSsheet("Sequencing Sheet")
    
    validateSampleSheetAgainstLIMS(lims, seq, limsMask, seqMask, projects="Maurano,CEGS")
