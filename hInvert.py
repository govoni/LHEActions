#!/usr/bin/python

import commands
import fileinput
import sys
import os

TESTING = 0
QUEUE = '8nm'


# copio nella cartella di arrivo meta' delle gen folder
# per la restante meta', chiamo la conversione

def runCommand (command, printIt = 0, doIt = 1) :
    if printIt : print ('> ' + command)
    if doIt : 
        commandOutput = commands.getstatusoutput (command)
        if printIt : print commandOutput[1]
        return commandOutput
    else :    print ('    command not run')
    return (1, 'command not run')
    

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


def prepareJob (tag, command) :
    filename = 'run_' + tag + '.job'
    f = open (filename, 'w')
    f.write ('cd /cvmfs/cms.cern.ch/slc6_amd64_gcc481/cms/cmssw/CMSSW_7_2_0\n')
    f.write ('eval `scram run -sh`\n')
    f.write ('cd -\n')
    f.write (command + '\n')
    f.close ()
    return filename


# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


def submitToLSF (tag, command) :
    print 'run : submitting jobs'
    jobFileName = prepareJob (tag, command)
    runCommand ('bsub -J ' + tag + ' -u pippopluto -q ' + QUEUE + ' < ' + jobFileName, 
                TESTING == 1, TESTING == 0)
    runCommand ('rm ' + jobFileName, TESTING == 1, TESTING == 0)


# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


if __name__ == "__main__":

    eoscmd = '/afs/cern.ch/project/eos/installation/cms/bin/eos.select' ;
    converter ='/afs/cern.ch/user/g/govoni/work/PHANTOM/phantom_at_cern/LHEActions/invertMuAndEleInLHE'

    inputFolder      = sys.argv[1]
    outputBaseFolder = os.path.abspath (sys.argv[2])
    outputUUFolder = outputBaseFolder + '/' + inputFolder + '_UE'
    outputEEFolder = outputBaseFolder + '/' + inputFolder + '_EU'
    inputFolder      = os.path.abspath (inputFolder)

    # find the gen folders
    res = runCommand ('ls ' + inputFolder)
    folders = [fol for fol in res[1].split () if 'gen' in fol and fol != 'gendir.scr']
    if (len (folders) != 100) :
        print 'PROBLEM: missing gen folders in ', inputFolder
        sys.exit (1)
    res = runCommand ('ls ' + outputBaseFolder)
    if (res[0] != 0) :
        print 'PROBLEM: missing output base folder ', outputBaseFolder
        sys.exit (1)
        
    res = runCommand ('ls ' + outputUUFolder, TESTING == 1, TESTING == 0)
    if (res[0] == 0) :
        print 'PROBLEM: output folder ', outputUUFolder, ' already existing'
        sys.exit (1)
    runCommand ('mkdir ' + outputUUFolder, TESTING == 1, TESTING == 0)
    
    res = runCommand ('ls ' + outputEEFolder, TESTING == 1, TESTING == 0)
    if (res[0] == 0) :
        print 'PROBLEM: output folder ', outputEEFolder, ' already existing'
        sys.exit (1)
    runCommand ('mkdir ' + outputEEFolder, TESTING == 1, TESTING == 0)
    
    # copy half of the events in the muons folder
    # ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
    
    runCommand ('cp ' + inputFolder + '/result ' + outputUUFolder, TESTING == 1, TESTING == 0)
    for i in range (len (folders) / 2) :
        runCommand ('cp -r ' + inputFolder + '/gen' + str (i + 1) 
                     + ' ' + outputUUFolder, TESTING == 1, TESTING == 0)

    # the electrons part
    # ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
    runCommand ('cp ' + inputFolder + '/result ' + outputEEFolder, TESTING == 1, TESTING == 0)
    for i in range (len (folders) / 2) :
        localInputFolder = inputFolder + '/gen' + str (i + 1 + len (folders) / 2)
        localOutputFolder = outputEEFolder + '/gen' + str (i + 1)
        runCommand ('mkdir ' + localOutputFolder,
                     TESTING == 1, TESTING == 0)
        submitToLSF (str (i+1),
                    converter + ' ' + localInputFolder + '/phamom.dat' 
                              + ' ' + localOutputFolder + '/phamom.dat')
    
