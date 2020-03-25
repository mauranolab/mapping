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


#getting warning /vol/mauranolab/mauram01/git/mapping/flowcells/lims.py:129: FutureWarning: `item` has been deprecated and will be removed in a future version
#https://stackoverflow.com/questions/57390363/pandas-item-has-been-deprecated
#TODO migrate to .values.item()?


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
        lims_gsheets_ID_file = open('/vol/mauranolab/flowcells/LIMS.gsheets_ID.txt', 'r')
        lims_gsheets_ID = lims_gsheets_ID_file.readlines()[0].rstrip("\n")
        lims_gsheets_ID_file.close()
        gc = pygsheets.authorize(service_file='/vol/mauranolab/flowcells/mauranolab.json')
        sh = gc.open_by_key(lims_gsheets_ID)
        wks = sh.worksheet_by_title(sheet)
        df = wks.get_as_df()
        df.fillna(value='', inplace=True) #Maybe pandas v1 problem? Or just change '' to None?
        #Remove per-FC headers and space between FCs (any row that is only empty lines)
        mask = ~df['Sample Name'].str.startswith('#') & (df!="").any(1)
    except Exception as e:
        #Doesn't print exception name right, e.g. if header is corrupted in google docs: "AttributeError: 'KeyError' object has no attribute 'argument'"
        print("WARNING could not load sheet " + sheet + " from google sheets: ", e.message, '\n', e.argument)
        wks = None
        df = None
    return wks, df, mask


#Verify consistency in common entries between Sample Sheet and LIMS sheets
#Projects argument suppresses less important inconsistencies unless project is on that comma-separated list
#TODO: ensure Trim empty except for Amplicon, Bait except for capture
def validateSampleSheetAgainstLIMS(lims, seq, limsMask, seqMask, projects="Maurano,CEGS"):
    projectList = projects.split(",")
    
    #TODO sort order
    print("#Check sequencing sheet for consistency with LIMS. Projects=", projects, sep="")
    print("Level", "Description", "Sample Name", "Sample #", "Key", "LIMS_value", "SampleSheet_value", sep=",")
    commonCols = set(lims.columns.values).intersection(set(seq.columns.values))
    numMissingSamples = 0
    numMultipleSamples = 0
    #Iterate through the sequencing sheet
    for seqRow in seq.index.values[seqMask]:
        curSeq = seq.iloc[seqRow]
        bs = curSeq['Sample #']
        SampleName = curSeq['Sample Name']
        curLims = lims[limsMask & lims['Sample #'].isin([bs])]
        numEntriesInLIMS = curLims.shape[0]
        if numEntriesInLIMS == 0:
            print("ERROR", "Can't find in LIMS!", SampleName, bs, "", "", "", sep=",")
            numMissingSamples += 1
        elif numEntriesInLIMS >= 2:
            print("ERROR", "found " + str(numEntriesInLIMS) + " entries in LIMS!", SampleName, bs, "", "", "", sep=",")
            numMultipleSamples += 1
        else:
            if projects=='' or curLims['Lab'].item() in projectList:
                for col in commonCols:
                    if curSeq[col] != curLims[col].item():
                        if curSeq[col] == "" or curLims[col].item() == "":
                            print("WARNING", "missing info", SampleName, bs, col, curLims[col].item(), curSeq[col], sep=",")
                        else:
                            print("ERROR", "inconsistent info", SampleName, bs, col, curLims[col].item(), curSeq[col], sep=",")
                for col in set(seq.columns.values):
                    if isinstance(curSeq[col], str) and curSeq[col] != curSeq[col].strip():
                        print("WARNING", "leading/trailing whitespace", SampleName, bs, col, "", curSeq[col], sep=",")
                for col in set(lims.columns.values):
                    if str(curLims[col].item()) != str(curLims[col].item()).strip():
                        print("WARNING", "leading/trailing whitespace", SampleName, bs, col, curLims[col].item(), "", sep=",")
    
    print(str(numMissingSamples) + " total samples missing from LIMS sheet")
    print(str(numMultipleSamples) + " total samples matching multiple entries in LIMS sheet")


#Pass in a data frame containing updates to make 
#Rows will be matched by 'Sample #', data in remaining columns will be transfered to the live copy of LIMS or Sample Sheet
#df must be the panda df mirror of wks, used to get coords efficiently
#TODO safer to make changes offline then sync? Not clear if this is implemented yet, unlink() seems to break get_value
#I think you do:
#wks.unlink()
#...
#wks.link(syncToCloud=True)
def updateSheetFromTable(wks, df, updates, commit=True):
    errors = 0
    colnames = df.columns.values
    excludedCols = ['Sample #'] #Don't update data in these columns
    for BS in updates['Sample #']:
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
                    newvalue = updates[updates['Sample #'] == BS][col].item()
                    if oldvalue != newvalue:
                        print(col + "=" + oldvalue + " => " + newvalue)
                    if commit:
                        wks.update_value(coords, newvalue)
    if errors > 0:
        print("ERROR: "+ errors + " found!")
