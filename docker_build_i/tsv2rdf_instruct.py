#coding: UTF-8
'''
Created on 2018/01/07
'''
import sys,os
from jinja2 import Template, Environment, FileSystemLoader
from optparse import OptionParser
import json
import numpy as np
import logging.handlers
import logging

RDF_TEMPLATE='templ_instruct.ttl'
EXEC_PATH = '/tmp'
MOUNT_PATH = '/mnt'
TEMPLATE_BODY = 'templ_instruct.ttl'
TEMPLATE_EVI = 'templ_instruct.ttl.evi'
TEMPLATE_PREFIX = 'templ_instruct.ttl.prefix'

class Message(object):
    def __init__(self, subject_no, data ):
        def get_db_name(publication):

            if publication.strip() == "Structural Genomics Consortium":
                db = "sgc"
            else:
                db = "pubmed"

            return db


        self.SNo = subject_no               #主語の連番
        self.ProtAUniprot = data["ProtA[uniprot]"]
        self.ProtBUniprot = data["ProtB[uniprot]"]
        self.ProtAOfficialSymbol = data["ProtA[Official Symbol]"]
        self.ProtBOfficialSymbol = data["ProtB[Official Symbol]"]
        self.PfamDomainA = data["Pfam-domainA"]
        self.PfamDomainB = data["Pfam-domainB"]
        self.DomainA = data["domainA"]
        self.DomainB = data["domainB"]
        self.DomainALoc = data["domainA_loc"].replace('"', '')
        self.DomainBLoc = data["domainB_loc"].replace('"', '')

        #Publications(PubMed or SGC)
        self.Publications = [ Publications(comma, publication, get_db_name(publication))
                             for comma, publication in enumerate(str(data["Publications"]).split(';')) ]   #PubMed

        #Publications(PDB)
        if data["PDB-structs"] is not np.nan:
            self.PDBs = [ Publications(1, publication, "pdb", row=data)
                             for evidence_no, publication in enumerate(data["PDB-structs"].split(';')) ]   #PDB
            self.Publications.extend(self.PDBs)
        else:
            #PDB-structsが空の場合がある
            pass

        return

class Publications(object):
    def __init__(self, comma, publication, db, row=None):

        self.Comma = comma      # 0: カンマなし >1:カンマあり
        self.DB = db            # "pdb" or "pubmed" or "sgc"
        if db == "sgc":
            self.ID = ""            #sgcの場合はIDなし
        else:
            self.ID = publication   #sgc以外の場合は文献のID

        if db == "pdb":
            if publication == row["Supporting_PDB_Structures"]:
                self.evidencelevel = "true"
            else:
                self.evidencelevel = "false"
        else:
            self.evidencelevel = ""

        return

def read_tsv(all, annotations, organism=None):
    '''
    @param all: Full structurally resolved interactome network with domain amino acid locations.
    @param annotations: Supplementary data.Domain-domain interactions shown to occur in specific protein-protein pairs within co-crystal structures. These annotations are combined from 3did and iPfam.
    '''
    import pandas as pd

    allDF=pd.read_csv(all, delimiter="\t")

    if annotations is not None:
        annotationsDF=pd.read_csv(annotations, delimiter="\t")
    else:
        annotationsDF=pd.DataFrame(columns=["Protein1","Protein2","Domain1","Domain2","Supporting_PDB_Structures"])

    #allのカラム名に合わせる
    annotationsDF.rename(columns={"Protein1" : "ProtA[uniprot]"}, inplace=True)
    annotationsDF.rename(columns={"Protein2" : "ProtB[uniprot]"}, inplace=True)
    annotationsDF.rename(columns={"Domain1" : "Pfam-domainA"}, inplace=True)
    annotationsDF.rename(columns={"Domain2" : "Pfam-domainB"}, inplace=True)

    mergedDF=pd.merge(allDF, annotationsDF, on=["ProtA[uniprot]", "ProtB[uniprot]","Pfam-domainA", "Pfam-domainB"], how="outer")

    #mergedDF.to_excel("debug_inst_%s.xlsx" % organism)

    #ProtA[Official Symbol]がNaN以外の行を抽出
    mergedDF = mergedDF[ ~mergedDF["ProtA[Official Symbol]"].isnull() ]

    return mergedDF


def main(fw, template, organism=None, all = None, annotations = None ):

    env = Environment(loader = FileSystemLoader(".", encoding='utf8'), autoescape = False)
    imageTemplate = env.get_template(template['body'])

    evidenceList=[]

    instruct_DF = read_tsv(all, annotations, organism=organism)

    for row in instruct_DF.iterrows():

        #Message Body
        msgObject = Message(row[0], row[1])
        evidenceList.extend([ '\t'.join([pub.DB, pub.ID, pub.evidencelevel]) for pub in msgObject.Publications])

        namespace = dict(message=msgObject)
        FeedContent = imageTemplate.render(namespace, organism=organism).encode('utf-8')

        fw.write(FeedContent)

    # テンプレートの読み込み(evidence)
    env = Environment(loader = FileSystemLoader(".", encoding='utf8'), autoescape = False)
    imageTemplate = env.get_template(template['evidence'])

    evidences_uq_wk = list(set(evidenceList))
    evidences_uq_wk2 = map(lambda x: x.split('\t'), evidences_uq_wk)
    evidences = [dict(DB=x[0],ID=x[1],evidencelevel=x[2]) for x in  evidences_uq_wk2]

    namespace = dict(Publications=evidences)
    FeedContent = imageTemplate.render(namespace).encode('utf-8')
    fw.write(FeedContent)

    return

if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option( "-c", type="string", help="config_file", dest="config", default= 'tsv2rdf_instruct.json')
    (options, args) = parser.parse_args()

    #ログ関連の設置
    argvs = sys.argv
    LOG_FILENAME = '%s.log' % argvs[0].split('.')[0]
    logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s',level=logging.DEBUG)
    logger = logging.getLogger()

    #ファイルログ
    file_log = logging.handlers.RotatingFileHandler(filename=LOG_FILENAME)
    file_log.setLevel(logging.DEBUG)
    file_log.setFormatter(logging.Formatter('%(asctime)s;%(levelname)s:%(message)s'))

    #ハンドラの設定
    logging.getLogger().addHandler(file_log)

    f = open(options.config ,"r")
    config = json.load(f)
    logger.info("Organism: %s" % ",".join([ k for k,v in config['organism'].items() ]))

    output_file = config['output_file']
    template = config['template']
    data_path = config['data_path']

    #出力先ファイル
    if os.path.exists(output_file) != True:
        output_file = os.path.join(MOUNT_PATH, output_file)
        fw = open(output_file, 'w')
    else:
        fw = open(output_file, 'w')

    #テンプレートファイルのチェック
    if os.path.exists(template['body']) != True:
        template['body'] = os.path.join(EXEC_PATH, template['body'])
        logger.info("Set template file path: %s" % template['body'])
    if os.path.exists(template['evidence']) != True:
        template['evidence'] = os.path.join(EXEC_PATH, template['evidence'])
        logger.info("Set template file path: %s" % template['evidence'])
    if os.path.exists(template['prefix']) != True:
        template['prefix'] = os.path.join(EXEC_PATH, template['prefix'])
        logger.info("Set template file path: %s" % template['prefix'])

    if os.path.exists(data_path) != True:
        data_path = os.path.join(MOUNT_PATH, data_path)
        logger.info("Set data_path: %s" % data_path)

    #prefixの出力
    fp = open(template['prefix'], 'r')
    prefixes = fp.readlines()
    for prefix in prefixes:
        fw.write(prefix)

    for k,v in config['organism'].items():

        if 'annotations' in v:
            annotations = os.path.join(data_path, v['annotations'])
        else:
            annotations = None

        logger.info("Processing organism: %s" % k)

        main(
            fw,
            template,
            organism = k,
            all = os.path.join(data_path, v['all']),
            annotations = annotations,
            )
    pass
