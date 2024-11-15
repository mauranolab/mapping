#!/bin/env python

#Functions for accessing Maurano Lab LIMS through pygsheets


#NB Need to get
#https://pygsheets.readthedocs.io/en/latest/authorization.html
#create project, service account key, share sheet with service account manually from google drive UI


import sys
import re
import argparse
import pandas as pd
import pygsheets
import os
import glob
from datetime import datetime

version="1.3"


#https://stackoverflow.com/questions/26987222/checking-whitespace-in-a-string-python/26987329
import string
def contains_whitespace(s):
    return True in [c in s for c in string.whitespace]


#https://stackoverflow.com/questions/16542074/what-is-the-inverse-of-date-toordinal-in-python
def parse_excel_date(exceldate):
    dt = datetime.fromordinal(datetime(1899, 12, 30).toordinal() + exceldate)
    return(dt.strftime("%Y-%m-%d"))


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
#wks is the live sheet, df is a panda df, and mask permits skipping rows as appropriate
def getLIMSsheet(sheet):
    try:
        #Read google sheets ID from the file system
        lims_gsheets_ID_file = open('/gpfs/data/isg_sequencing/LIMS.gsheets_ID.txt', 'r')
        lims_gsheets_ID = lims_gsheets_ID_file.readlines()[0].rstrip("\n")
        lims_gsheets_ID_file.close()
        
        gc = pygsheets.authorize(service_file='/gpfs/data/isg_sequencing/mauranolab.json')
        sh = gc.open_by_key(lims_gsheets_ID)
        wks = sh.worksheet_by_title(sheet)
        #get_as_df() doesn't permit setting date_time_render_option separately, so value_render="UNFORMATTED_VALUE" seems to result in dates being output as serials. UNFORMATTED_VALUE is needed to avoid thousand separators, and decimal formatting inconsistencies.
        df = wks.get_as_df(value_render="UNFORMATTED_VALUE")
        df.fillna(value='', inplace=True) #Maybe pandas v1 problem? Or just change '' to None?
        #Mask per-FC headers and space between FCs (any row that is only empty lines)
        mask = ~df['Sample Name'].astype(str).str.startswith('#') & (df!="").any(axis=1)
    except Exception as e:
        #Doesn't print exception name right, e.g. if header is corrupted in google docs: "AttributeError: 'KeyError' object has no attribute 'argument'"
        print("WARNING could not load sheet " + sheet + " from google sheets: ", e.message, '\n', e.argument)
        wks = None
        df = None
    return wks, df, mask


#Verify consistency in common entries between Sample Sheet and LIMS sheets
#Projects argument suppresses less important inconsistencies unless project is on that comma-separated list
def validateSampleSheetAndLIMS(lims, seq, limsMask, seqMask, projects="", runnames="", sampleids="", madeby=""):
    projectList = set(filter(None, projects.split(",")))
    sampleidList = set(filter(None, sampleids.split(",")))
    madebyList = set(filter(None, madeby.split(",")))
    #Implement runname filtering by making a pass through the sequencing sheet and adding BS numbers from the specified run to sampleidList
    if runnames != "":
        runnameList = runnames.split(",")
        
        curRunname = None
        #Mask space between FCs (any row that is only empty lines)
        for seqRow in seq.index.values[(seq!="").any(axis=1)]:
            curSeq = seq.iloc[seqRow]
            if curSeq['Sample Name'] == "#Run name":
                curRunname = curSeq['Sample #']
            elif curRunname in runnameList:
                #Translate runnameList into sampleidList
                #BUGBUG does union rather than intersect
                #Skip remaining per-FC headers with leading "#"
                if re.match("^#", str(curSeq['Sample Name'])) is None:
                    sampleidList.add(curSeq['Sample #'])
    
    
    print("Level", "Description", "Sample Name", "Sample #", "Key", "LIMS_value", "SampleSheet_value", sep="\t")
    
    print("#Check sequencing sheet for consistency with LIMS. Projects=", projects, sep="", file=sys.stderr)
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
            if (len(projectList)==0 or curLims['Lab'].values.item() in projectList) and (len(sampleidList)==0 or bs in sampleidList) and (len(madebyList)==0 or curLims['Made By'].values.item() in madebyList):
                for col in commonCols:
                    if curSeq[col] != curLims[col].values.item():
                        #BUGBUG I tried switching to "".__eq__(myString) but somehow it always evaluates to true. https://stackoverflow.com/questions/9573244/how-to-check-if-the-string-is-empty#
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
    
    print("#Check LIMS for integrity. Projects=", projects, sep="", file=sys.stderr)
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
        
        #Only check certain metadata for specified projects to avoid excess verbiage
        if (len(projectList)==0 or curLims['Lab'] in projectList) and (len(sampleidList)==0 or bs in sampleidList) and (len(madebyList)==0 or curLims['Made By'] in madebyList):
            for col in ["Sample Name"]:
                if re.match("[\+%\(\)\"\'\/\. ]", str(curLims[col])) is not None:
                    print("WARNING", "illegal characters in LIMS", SampleName, bs, col, curLims[col], "", sep="\t")
                #Hyphens are only allowed for ChIP-seq samples, where 1 is required
                samplename_hyphen_re = re.findall("\-", str(curLims[col]))
                if str(curLims['Sample Type']) in ['ChIP-seq', 'CUT&RUN']:
                    if len(samplename_hyphen_re) !=1:
                        print("WARNING", "missing ChIP-seq epitope", SampleName, bs, col, curLims[col], "", sep="\t")
                else:
                    if len(samplename_hyphen_re) !=0:
                        print("WARNING", "illegal characters in LIMS", SampleName, bs, col, curLims[col], "", sep="\t")
            
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
                    print("ERROR", "duplicate custom reference", SampleName, bs, "Genetic Modification/Custom Reference", sorted(geneticModifications & customReferences), "", sep="\t")
            
            requiredColsBySampleType = { "DNA Capture": ["Parent Library", "Bait set", "Pool ID"] }
            for sampleType in requiredColsBySampleType:
                for col in requiredColsBySampleType[sampleType]:
                    if curLims["Sample Type"] in [sampleType]:
                        if curLims[col] == "":
                            print("WARNING", col + " column is required for sample type " + sampleType, SampleName, bs, col, curLims[col], "", sep="\t")
            
            colsSpecificToSampleType = { "R1 Trim (P5)": ["Amplicon", "Transposon DNA", "Transposon RNA", "Transposon iPCR", "Transposon 10xRNA"], "R2 Trim (P7)": ["Amplicon", "Transposon DNA", "Transposon RNA", "Transposon iPCR", "Transposon 10xRNA"], "Bait set" : ["DNA Capture", "Amplicon", "Pool", "Fiber-seq"] }
            for col in colsSpecificToSampleType:
                if curLims[col] != "":
                    if curLims["Sample Type"] not in colsSpecificToSampleType[col]:
                        print("WARNING", col + " column is not allowed for sample type " + curLims["Sample Type"], SampleName, bs, col, curLims[col], "", sep="\t")
    
    
    print("", file=sys.stderr)
    print("#Clone info", file=sys.stderr)
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
#The paranoid option verifies the BS and value to update directly in wks rather than trusting df, mainly to guard against coordinate arithmetic errors, etc
#TODO safer to make changes offline then sync? Looks like link()/unlink() has been deprecated
def updateSheetFromTable(wks, df, updates, commit=True, paranoid=True):
    if any(updates['Sample #'].duplicated()):
        print("ERROR: found duplicate Sample #s in update table!")
        return
    errors = 0
    colnames = df.columns.values
    excludedCols = ['Sample #'] #Don't update data in these columns
    for BS in updates['Sample #']:
        if re.match("^BS[0-9]+[A-Z]$", BS) is None:
            print("WARNING: empty or invalid BS provided '" + BS + "'")
            continue
        rows = df.index[df['Sample #'] == BS].values
        #I think row + 2 because of 0 -> 1 based indexing plus the header row
        for row in rows + 2:
            row = int(row) #Needed to explicitly cast to int after 2024 bigpurple migration
            print("\nRow " + str(row) + ":" + BS)
            for col in set(colnames).intersection(set(updates.columns.values)).difference(excludedCols):
                #Need to translate row back to prior coordinates
                oldvalue = df.iloc[row-2][col]
                if col in ['QC Gel Date', 'Bioanalyzer Run Date', 'Date Library Finished']:
                    if oldvalue != "":
                        oldvalue = parse_excel_date(oldvalue)
                coords = (row, df.columns.get_loc(col)+1)
                newvalue = updates[updates['Sample #'] == BS][col].values.item()
                if oldvalue != newvalue:
                    if paranoid:
                        wksBS = str(wks.get_value((row, df.columns.get_loc('Sample #')+1)))
                        if BS != wksBS:
                            errors += 1
                            print("ERROR: " + wksBS + " does not match expected sample " + BS)
                        
                        wksoldvalue = wks.get_value(coords)
                        #BUGBUG seems to fail incorrectly when changing int fields with existing data? Also confused by date format mismatch.
                        if oldvalue != wksoldvalue:
                            errors += 1
                            print("ERROR: " + BS + " wks value " + str(wksoldvalue) + " does not match df value " + str(oldvalue))
                    
                    print(col + "=" + str(oldvalue) + " => " + str(newvalue))
                    if commit:
                        wks.update_value(coords, newvalue)
    if errors > 0:
        print("ERROR: "+ str(errors) + " found!")


#Identify FCs that contain a given list of samples
#NB does not actually check whether changes were made
#TODO could be done with less validation using indices
def findFCsForSamples(seq, updates):
    affectedFCs = []
    curFC = None
    numEmptyLines = 0
    for seqRow in seq.index.values:
        curSeq = seq.iloc[seqRow]
        if (curSeq == "").all():
            numEmptyLines += 1
            #clear curFC if we pass two newlines
            if numEmptyLines >= 2:
                curFC = None
            continue
        else:
            numEmptyLines = 0
            if curSeq['Sample Name'] == "#Barcode":
                curFC = curSeq['Sample #']
            elif re.match("^#", str(curSeq['Sample Name'])) is None:
#                print(curSeq['Sample #'], curFC)
                if curSeq['Sample #'] in updates['Sample #'].values:
                    if curFC is not None:
                        affectedFCs.append(curFC)
#                    print("Found:", curSeq['Sample #'], "on FC", curFC)
    print(set(affectedFCs), sep=",")


#Search and replace functionality for FC info
#BUGBUG I don't think this is safe to run twice in a given session
def replaceFCinfo(seqWks, seq, key, value, newvalue, commit=False):
    fcMask = seq['Sample Name'].astype(str).str.startswith('#')
    curFC = None
    for seqRow in seq.index.values[fcMask]:
        curSeq = seq.iloc[seqRow]
        if curSeq['Sample Name'] == "#Barcode":
            #BUGBUG FC is customarily after run name, so it won't be set for the first line
            curFC = curSeq['Sample #']
        elif curSeq['Sample Name'] == "#"+key:
            #I think row + 2 because of 0 -> 1 based indexing plus the header row
            coords = (seqRow+2, seq.columns.get_loc('Sample #')+1)
            oldvalue = seqWks.get_value(coords)
            if (value is None or oldvalue == value) and oldvalue!=newvalue:
                print("Row ", seqRow+2, " (", curFC, "): ", key, "=", str(oldvalue), " => ", str(newvalue), sep="")
                if commit:
                    seqWks.update_value(coords, newvalue)


#Remove empty FC info records matching a given key
def removeEmptyFCinfo(seqWks, seq, key, commit=True):
    fcMask = seq['Sample Name'].astype(str).str.startswith('#')
    #Iterate in reverse order since updating row index to account for deleted rows does not seem to work
    for seqRow in reversed(seq.index.values[fcMask]):
        curSeq = seq.iloc[seqRow]
        if curSeq['Sample Name'] == "#"+key and "".__eq__(curSeq['Sample #']):
            print("Removing row ", seqRow, ": ", curSeq['Sample Name'], sep="")
            if commit:
                #I think row + 2 because of 0 -> 1 based indexing plus the header row
                seqWks.delete_rows(int(seqRow+2), 1)


###Command line operation
if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog = "lims", description = "Utilities for Maurano Lab LIMS implemented in google sheets. Validates or updates LIMS entries.", allow_abbrev=False)
    
    validate_parser = parser.add_argument_group()
    validate_parser.add_argument("--validate", action = "store_true", default = False, help = "Run LIMS and Sample Sheet validation [%(default)s]")
    validate_parser.add_argument("--runnames", action = "store", type = str, default="", help = "Only process samples contained in these run names (multiple run names separated by comma, done as string matching)  [%(default)s]")
    validate_parser.add_argument("--sampleids", action = "store", type = str, default="", help = "Only process samples matching this BS number (multiple samples separated by comma, done as string matching) [%(default)s]")
    validate_parser.add_argument("--projects", action = "store", type = str, default="Maurano,CEGS", help = "Only process samples in these projects (multiple projects separated by comma) [%(default)s]")
    validate_parser.add_argument("--madeby", action = "store", type = str, default="", help = "Only process samples made by these people (multiple names separated by comma) [%(default)s]")
    
    update_parser = parser.add_argument_group()
    update_parser.add_argument("--update", action = "store", default = None, type=str, help = "Update LIMS and Sample Sheet based on tab-delimited input file [%(default)s]")
    update_parser.add_argument("--nocommit", action = "store_true", default = False, help = "Perform update but do not commit any differences [%(default)s]")
    
    parser.add_argument('--version', action='version', version='%(prog)s ' + version)
    
    #argparse does not set an exit code upon this error
    #https://stackoverflow.com/questions/5943249/python-argparse-and-controlling-overriding-the-exit-status-code
    try:
        args = parser.parse_args()
    except argparse.ArgumentError as exc:
        print(exc.message, '\n', exc.argument)
        sys.exit(2)
    
    #Load google sheets
    limsWks, lims, limsMask = getLIMSsheet("LIMS")
    seqWks, seq, seqMask = getLIMSsheet("Sequencing Sheet")
    
    #Perform work as directed on command line
    if args.validate:
        validateSampleSheetAndLIMS(lims, seq, limsMask, seqMask, projects=args.projects, runnames=args.runnames, sampleids=args.sampleids, madeby=args.madeby)
    
    if args.update is not None:
        if args.nocommit:
            print("Trial update only -- no edits will be committed.\n")
        
        updates = pd.read_csv(args.update, na_filter=False, sep="\t")
        print("Updating LIMS and Sequencing Sheet based on:")
        print(updates)
        
        print()
        print("Updating LIMS:")
        updateSheetFromTable(limsWks, lims, updates, commit=not args.nocommit)
        
        print()
        print("Updating Sequencing Sheet:")
        updateSheetFromTable(seqWks, seq, updates, commit=not args.nocommit)
        
        print()
        print("FCs containing the samples in this update")
        findFCsForSamples(seq, updates)

