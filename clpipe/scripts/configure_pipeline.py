import sys
import os
import shutil
import argparse

import ruamel.yaml
yaml_ruamel = ruamel.yaml.YAML()
yaml_ruamel.indent(mapping=4, sequence=4, offset=2)

import json


# --------------------------------------------------------------------
# Functions used to read access specific data in a yaml file

def GetYamlDocumentSection(sectionID, filename=None):

    # Get configuration and setup all values
    if filename==None: filename = input_yaml_file
    with open(filename) as file:
        raw_text = file.read()
    generic_config_docs =  yaml_ruamel.load_all(raw_text)

    for doc in generic_config_docs:
        if doc and "id" in doc and doc["id"]==sectionID: return doc

    return None

def GetYamlDocumentSectionList(filename=None):

    if filename==None: filename = input_yaml_file

    idList=[]
    with open(filename) as file:
        raw_text = file.read()
    generic_config_docs =  yaml_ruamel.load_all(raw_text)

    for doc in generic_config_docs:
        if doc and "id" in doc : idList.append(doc["id"])

    return idList


# --------------------------------------------------------------------
# Usefull functions 

def ExtractSortedConfigureValuesFromDict(valueDict):
    
    keyValues = list(valueDict)
    keyValues.sort(key=len)
    keyValues.reverse()
    #print(keyValues)

    return keyValues


# --------------------------------------------------------------------
# Read and initalize the configuration paramaters from the configure section

def ReplaceKeyValuesByRealValues(keyList,valueDict):

    # lopp over keyValues
    bInitialization = True
    while(bInitialization):
        bInitialization = False
        for key in keyList:
            for k,v in valueDict.items():
                if key in v:
                    v_new=v.replace(key,valueDict[key])
                    valueDict[k]=v_new

    return

def ReadAndSetupConfigurationSection(overwriteValues=None):

    # Read configuration values from configuration file
    dataDict = GetYamlDocumentSection("configure")

    # Overwrite values if some are defined from python script input parameters
    if isinstance(overwriteValues,dict):
        for k,v in overwriteValues.items():
            oldValue = None
            if k in dataDict: oldValue = dataDict[k]
            print(f"OVERWRITE config : {k} {v} vs {oldValue}")
            dataDict[k]=v

    if "site" in dataDict:
        siteName = dataDict["site"]
    else:
        print("No site name defined in id:confirgure section")
        sys.exit()

    # Get list of parameters defined in id:configure
    configKeys = list(dataDict)
    configKeys.remove("id")
    if "name" in configKeys: dataDict["name"] = dataDict["name"].replace(" ","_")

    # Define config parameter values
    keyValueDict={}
    for key in configKeys:
        v = dataDict[key]
        if isinstance(v,dict) :
            if not siteName in v or v[siteName]==None: 
                raise Exception(f"Configuration : {siteName} value not defined for {key}")
            v = v[siteName]
        keyValueDict["$"+key.upper()] = v

    # add local env parameters
    keyValueDict["$PWD"] = os.environ["PWD"]
    keyValueDict["$HOME"] = os.environ["HOME"]

    # Extract config parameter keys (sorted by string length)
    keyValues = ExtractSortedConfigureValuesFromDict(keyValueDict)

    # Initialize config parameters - replace $VALUES by real values
    ReplaceKeyValuesByRealValues(keyValues, keyValueDict)

    # Final configuration parameters
    print("---- Configuration parameters :")
    import pprint
    pprint.pprint(keyValueDict)
    print("-"*25)

    return keyValueDict


# --------------------------------------------------------------------
# Read and initalize the configuration paramaters from the configure section

def CreateConfigurationFilesDirectory(pDict):
    """ Create the directory to store the configured files """

    configDir = pDict["$NAME"]
    if not os.path.isdir(configDir): os.makedirs(configDir)

    return configDir

def ConfigureFileBasedOnParamDict( filename, paramDict):
    """ Replace the configuration values by their real values """

    # get the dict keys sorted by their lengths
    keyList = ExtractSortedConfigureValuesFromDict(paramDict)

    # Use linux sed command to replace config values by real values
    for k in keyList:
        v=paramDict[k]
        v_new=v.replace("/","\/")
        v_new=v_new.replace("\n","\\n")
        cmd=f"sed -i 's/{k}/{v_new}/g' {filename}"
        #print("CMD "+cmd)
        os.system(cmd)
    
    return


def CreateAndConfigureYamlConfigFile(pDict):
    """ Create and configure the yaml config file """

    ymlConfigFileGen = pDict["$CONFIG_FILE_GENERIC"]
    ymlConfigFile = pDict["$COMPUTATIONDIR"]+"/"+pDict["$CONFIG_FILE"]

    shutil.copy(ymlConfigFileGen,ymlConfigFile)
    ConfigureFileBasedOnParamDict( ymlConfigFile, pDict)

    return

def CreateJobConfigurationFiles(paramDict):
    """ Split the concatenated job file, create and configure the TXpipe, FireCrown, etc files  """

    # Get yaml job file
    jobfilename = None
    if "$JOBCONCATYAMLFILE" in paramDict: jobfilename=paramDict["$JOBCONCATYAMLFILE"]

    # Get id list from yaml concatenated file
    idList = GetYamlDocumentSectionList(jobfilename)

    # Create working directory
    dirName = CreateConfigurationFilesDirectory(paramDict)

    # Create and configure yaml files for each job
    for jobId in idList:
        if jobId=="configure": continue

        # Create simple yaml file
        data = GetYamlDocumentSection(jobId, jobfilename)
        filename = f"{dirName}/{jobId}.yml"
        yaml_ruamel.dump(data,open(filename, 'w'))

        # Setup file values based on configuration parameters
        ConfigureFileBasedOnParamDict(filename,paramDict)

    # Create yaml config file  
    CreateAndConfigureYamlConfigFile(paramDict)

    return


def CreateBatchSubmissionScripts(paramDict):
    """ Create and configure the scripts used to submit the jobs """

    # Retrieve batch setup from configration section and batch #SBATCH param
    if "$BATCH_ID" in paramDict: 
        batchId = paramDict["$BATCH_ID"]
    else:
        return
    print(f"BATCH {batchId}")

    # Read batch configutaion parameters from batch_param.yaml file
    filename = "./batch_param.yml"
    if "$BATCH_PARAM_FILE" in paramDict:
        filename = paramDict["$BATCH_PARAM_FILE"]
    batchParam = GetYamlDocumentSection(batchId, filename)
    print(batchParam)
    
    # Get site value
    sitename = paramDict["$SITE"]

    for k,v in batchParam.items():
        if k=="id": continue
        if sitename in v and v[sitename]!=None: 
            param = v[sitename]
            paramStr = "\n".join([f"#SBATCH {x.strip()}" for x in param.split(";") if len(x)>0])
            paramDict[f"${k.upper()}"]=paramStr

    # Create batch files
    tmp = paramDict["$BATCH_FILES"]
    batchfiles=[x for x in tmp.split(";")]
    for filename in batchfiles:
        filename_res = paramDict["$COMPUTATIONDIR"]+"/"+filename.replace("generic",sitename).split("/")[-1]
        shutil.copy(filename, filename_res)
        ConfigureFileBasedOnParamDict(filename_res, paramDict)

    return



if __name__=="__main__":

    # Input parameters : configuration file and parameter dict
    input_yaml_file=sys.argv[1]

    overwriteDict={}
    if len(sys.argv)==3:
        overwriteDict = json.loads(sys.argv[2])

    # parser = argparse.ArgumentParser()
    # parser.add_argument('-f', '--file', type=str)
    # parser.add_argument('-d', '--mydict', type=json.loads)
    # args = parser.parse_args()
    # overwriteDict = args.mydict  # This will now return a dictionary

    # Read and initialize id:configure section
    configValuesDict = ReadAndSetupConfigurationSection(overwriteDict)

    # Create job configurations files - split file into small yamlfile + config files
    CreateJobConfigurationFiles(configValuesDict)

    # Build batch submission scripts
    CreateBatchSubmissionScripts(configValuesDict)

    # Build parsl workflow


