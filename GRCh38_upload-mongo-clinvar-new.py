#!/usr/bin/python
# -*- coding: utf-8 -*-
#================================================================================#
# Python script to import a left alligned VCF format file into mongodb.          #
# Assumes a FEATURE field as the first field in INFO in JSON format.             #
#================================================================================#
import pymongo
import sys
import util
from pymongo import Connection
from pymongo.collection import Collection
import re

#============================================================#
#                 initialize globals                         #
#============================================================#
#----- script arguments -----#
args = None

#----- VCF CLN delimeters -----#
DELIM_COMMA  = ","
DELIM_BAR    = "|"
DELIM_EQUAL  = "="

#----- VCF CLN field names -----#
CLNACC    = "CLNACC"
CLNALLE   = "CLNALLE"
CLNDBN    = "CLNDBN"
CLNDSDB   = "CLNDSDB"
CLNDSDBID = "CLNDSDBID"
CLNHGVS   = "CLNHGVS"
CLNORIGIN = "CLNORIGIN"
CLNSIG    = "CLNSIG"
CLNSRC    = "CLNSRC"
CLNSRCID  = "CLNSRCID"

#----- list of clinical keys to process -----#
# these are the fields that can have multiple values.
ClinicalKeys = [ CLNACC, CLNALLE, CLNDBN, CLNDSDB, CLNDSDBID, CLNHGVS, CLNORIGIN, CLNSIG, CLNSRC, CLNSRCID ]

#----- flags for CLNALLE field -----#
CLNALLE_NO_ALLELE_WAS_FOUND = "-1"

#----- flags for no allele -----#
NO_ALTERNATIVE_ALLELE = "A"

#----- mongo database name for annotations -----#
DATABASE_NAME = "AnnotationSource"

#----- create map from clinical significance (CLNSIG) field from a number to string -----#
clnsigMap = { '0':   'Uncertain significance',
              '1':   'not provided',
              '2':   'Benign',
              '3':   'Likely benign',
              '4':   'Likely pathogenic',
              '5':   'Pathogenic',
              '6':   'drug response',
              '7':   'histocompatibility',
              '255': 'other'
           }

#============================================================#
#                       parse_args                           #
#============================================================#
# parse the input arguments.
# NOTE: this is not currently used.

def parse_args():
    global args
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf-file",            help="name of the VCF file")
    parser.add_argument("--collection-name",     help="name of the databse collection to create")
    parser.add_argument("--collection-version",  help="name of the databse collection version")

    #----- check arguments -----#
    if len(sys.argv) < 6: 
        print """
        Usage: GRCh38_upload_mongo_clinvar-new.py -m <mapping file>  <vcffile1>,<vcffile2>... <dbcollectionname> <dbcollectionversion> 
	Example: python GRCh38_upload_mongo_clinvar-new.py -m GCF_000001405.28.assembly.txt clinvar_20160203.vcf,clinvar_20160203_papu.vcf GRCh38_clinvar 20160203
        """
        sys.exit()

#============================================================#
#                  parseVariantFields                        #
#============================================================#
# parse the VCF record fields. This will consist of information
# fields and clinical fields CLNSIG, CLNACC, CLNDBN, CLNDSDBID 
# and CLNDSDB. these fields may have multiple values depending 
# on if clinical alternative alleles are present and/or multiple 
# assertions of clinical significance as refenced in 
# ClinVar (http://www.ncbi.nlm.nih.gov/clinvar/):
#
#    1) CLNALLE=1     CLNSIG=255
#    2) CLNALLE=1     CLNSIG=5|255   multiple assertions of clinical significance in ClinVar
#    3) CLNALLE=1,2   CLNSIG=5,5     multiple alleles
#    4) CLNALLE=1,2   CLNSIG=5,5|5
#
# The non-clinical fields are returned in a single dictionary while
# clinical fields with multiple values are retured as a list of dictionarires.

def parseVariantFields( observed_alleles, info, filterkey ): 

    variant_fields = info.split(';')
    variant_data = {}
    clinical_data_list = []
    clinical_alleles = []

    #----- parse variant fields -----#
    for field in variant_fields:
        if ( DELIM_EQUAL not in field ):
            continue
        key = field.split(DELIM_EQUAL)[0]
        value = field.split(DELIM_EQUAL)[1]

        if ( filterkey == None ):
            filterkey = key;

        # create a dictionary for each clinical allele
        if ( key == CLNALLE ):
           clinical_alleles = value.split(DELIM_COMMA)
           n = len(clinical_alleles)
           for i in xrange(n):
              clinical_data_list.append( {} ) 

        if ( key in ClinicalKeys ):
           addKeyValue( clinical_data_list, key, value )
        else:
           variant_data[key] = value
           
    #__for field in variant_fields__ 

    return variant_data, clinical_data_list, clinical_alleles, filterkey 

#============================================================#
#                  createVariantRecord                       #
#============================================================#
# create a variant record for entry into the database. this
# will create a list of dictionaries for clinical data separated
# by '|'s. 

def createVariantRecord( id, reference_allele, variant_data, clincal_allele, clinical_data, observed_alles ):

    #----- set obsevered allele -----#
    if ( clincal_allele == CLNALLE_NO_ALLELE_WAS_FOUND ):
       observed_allele = NO_ALTERNATIVE_ALLELE 
    elif ( clincal_allele == 0 ):
       observed_allele = reference_allele 
    else:
       i = int(clincal_allele) - 1
       observed_allele = observed_alles[i]

    records = []
    num_fields = 1

    #----- check for multiple CLNACC fields -----#
    if ( CLNACC in clinical_data ):
       clnacc_data = clinical_data[CLNACC]

       if ( not isinstance(clnacc_data, basestring) ):
          num_fields = len(clnacc_data)

    #----- add clinical fields -----#
    for i in xrange(num_fields):
       record = variant_data.copy()
       record["id"] = id
       record["r"] = reference_allele
       record["o"] = [ observed_allele ]

       for key in clinical_data:
          if ( key == CLNALLE ):
             record[key] = clincal_allele 
          else:
             data = clinical_data[key]
             if ( not isinstance(data, basestring) and (len(data) == num_fields) ):
                record[key] = data[i]
             else:
                record[key] = data

       records.append( record ) 

    return records

#============================================================#
#                      addKeyValue                           #
#============================================================#
# add a key-value pair for clinical data to the list of records. 
# 'value' may contain several values separated by commas or |'s. 
# if key = CLSIG then map its numeric value to a descriptive
# string.

def addKeyValue( data_list, key, value ):

    if ( key == CLNDBN ):
       value = value.replace("\\x2c","")

    if ( DELIM_COMMA in value ):
       value_list = value.split( DELIM_COMMA )
    else:
       value_list = [ value ] 

    for i in xrange(len(value_list)):
       v = value_list[i]
       data = data_list[i]

       if ( DELIM_BAR in v ):
          listv = v.split(DELIM_BAR)
          if ( key == CLNSIG ):
             sig_list = []
             for sig in listv:
                sig_list.append( convClinSig(sig) ) 
             data[key] = sig_list
          elif (key == CLNACC ):
             acc_list = []
             for acc in listv: 
                 splitACC=acc.split(".")[0]
                 acc_list.append(splitACC)
                 #print acc_list
             data[key] = acc_list
          else:
             data[key] = listv
       else:
          if ( key == CLNSIG ):
             data[key] = convClinSig( v )
          elif ( key == CLNACC ):
             data[key] = v.split(".")[0]
          else:
             data[key] = v

#============================================================#
#                        convClinSig                         #
#============================================================#
# convert a numerical clinical significance value into a 
# descriptive string.

def convClinSig( value ):
    if ( value not in clnsigMap ):
       sys.stderr.write( "**** ERROR: Unknown significance value \"%s\" \n" % value ) 
       return clnsigMap["255"]
    else:
       return clnsigMap[value]

#============================================================#
#                        main                                #
#============================================================#
# connect to mongo database, parse VCF records and add them
# to mongo.

def main():

    #----- process arguments -----#
    args = util.getArgs(sys.argv)
    mapping = util.readmapping(args["mapfile"])

    #----- connect to the mongo database -----#
    c = Connection(args["host"], args["ip"])
    db = c[DATABASE_NAME]
    db.drop_collection(args["table_name"])
    vcf = db[args["table_name"]]
    sys.stdout.write( "Created mongoDB table %s \n" % args["table_name"] )
    debug = False
    #debug = True

    #----- read vcf file and process variant records -----#
    count = 0
    filterkey = None
    hg19_filter_count = 0
    files = args["datafile"].split(",")
    for ff in files:
        f = open(ff, "r")

        for line in f:
            line = line.strip()

            if line[0] == '#':
                continue

            fields = line.split("\t")

            if len(fields) < 8:
                sys.stdout.write("%s\n" % line )
                sys.stdout.write("Invalid vcf, Not enough columns in the provided vcf file. \n" )
                sys.exit(1)
        
            #----- extract basic info: chromosome id, position, etc. -----#
            tchr = util.getValue(0, fields)
            if tchr in mapping:
                tchr = mapping[tchr]
                if tchr == "na":
                    print "Unrecognized chromosome: " + line
                    continue
            chr = util.getChr(tchr)
            ncbi = util.getNcbi(tchr)
            ct = util.getChrType(tchr)
            pos = int(util.getValue(1,fields))
            id = util.getValue(2,fields)
            reference_allele = util.getValue(3,fields)
            observed_alleles = []
            listv = []
            observed_alleles = util.getValue(4,fields).split(",")
            info = util.getValue(7,fields)

            #---------- parse informaion fields ----------#
            # returns list of possibly multiple dicts for 
            # clinical fields and a single dict for the rest.
            variant_data, clinical_data_list, clinical_alleles, filterkey = parseVariantFields( observed_alleles, info, filterkey )

            #---------- parse genotype field ----------#
            if len(fields) == 10:
                if "GT" in fields[8]:
                    GTFormat = util.getValue(8,fields)
                    GTValue = util.getValue(9,fields)
                    GT_record = {}
                    k = GTFormat.split(':')[0] 
                    v = GTValue.split(':')[0] 
                    variant_record[k] = v;

            #----- add record to database -----#
            if ( debug ): 
                sys.stdout.write("########## chr=%s pos=%s obsA=%s ######### \n" % (chr, pos, observed_alleles) )
                sys.stdout.write(">>> num_clin_alleles=%s \n" % len(clinical_alleles) )
            database_entry = { "_id" : { "c" : chr, "p" : pos } }
            if ncbi:
                database_entry["_id"]["ncbi"] = ncbi
            if ct:
                database_entry["_id"]["ct"] = ct
            count = count + 1

            for i in xrange(len(clinical_alleles)):
                clincal_allele = clinical_alleles[i]
                clinical_data = clinical_data_list[i]
                variant_records = createVariantRecord( id, reference_allele, variant_data, clincal_allele, clinical_data, observed_alleles )
                # add multiple features 
                for variant_record in variant_records:
                    pushUnique = { "$addToSet":{ "f" : variant_record } }
                    vcf.update( database_entry, pushUnique, True, False)
                    if ( debug ): 
                       sys.stdout.write("\n")
                       sys.stdout.write(">>> variant_record=%s \n" % (variant_record) )
                #__for variant_record in variant_records__
            #__for i in xrange(len(clinical_alleles))__
        #__for line in f__

    #----- add a meta data entry -----#""
    sys.stdout.write("Adding meta record \n")
    meta = { "_id" : "meta",
             "level" : "variant",
             "name" : "clinvar",
             "type" : "CLINVAR",
             "version" : 20160203,
             "links" : ["http://www.ncbi.nlm.nih.gov/clinvar/{CLNACC}"],
             "reference" : ["GRCh38"],
             "desc" : "ClinVar",
             "buildVer" : ["5.2"]
           }

    if filterkey == None:
        filterkey = "CLNSIG"
    meta["filterKeys"] = [filterkey]

    vcf.insert(meta)
    sys.stdout.write("Inserted meta record \n")

    vcf.ensure_index([ ("_id.c", pymongo.ASCENDING),
                       ("_id.ncbi", pymongo.ASCENDING),
                       ("_id.ct", pymongo.ASCENDING),
                       ("_id.p", pymongo.DESCENDING)
                     ])  
    version = [] 
    ins_count = vcf.count() - 1
    sys.stdout.write(" ********** Processed %d records in file ********** \n" % count )
    sys.stdout.write(" ********** Inserted %d records in table ********** \n" % ins_count )
    sys.stdout.write(" ********** Filtered for missing hg19 coords %d records in table ********** \n" % hg19_filter_count ) 

#============================================================#
#                    main program                            #
#============================================================#

if __name__ == '__main__':
    main()

