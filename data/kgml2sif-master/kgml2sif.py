#!/usr/bin/env python3
#### LICENSE ###################################################################
# Copyright (C) 2019-2020 Arnaud Poret
# This work is licensed under the GNU General Public License.
# To view a copy of this license, visit https://www.gnu.org/licenses/gpl.html
#### IMPORTS ###################################################################
import argparse
import copy
import os.path
import re
import xml.etree.ElementTree
#### PARSER ####################################################################
parser=argparse.ArgumentParser(
description="""
Convert KGML encoded KEGG human pathways to the SIF file format.
""",
epilog="""
Currently:
    * the following node types are considered:
        * gene (in KGML, genes also stand for the corresponding proteins/gene
          products)
        * compound
        * group (in KGML, groups stand for complexes)
    * the following relation types are considered (in KGML, relations stand for
      edges):
        * PPrel: protein-protein relations
        * GErel: protein-gene relations (i.e. gene expression regulations)
        * PCrel: protein-compound relations
        * ECrel: enzyme-enzyme relations sharing a common compound (as product
          for the first enzyme and as substrate for the second one)

In the output SIF file:
    * edges for which the name is missing are automatically named "unknown"
    * multiedges, if any, are merged if possible, otherwise they are left inside
      the output SIF file and kgml2sif warns about them
    * the non KGML "membership" edges are added in order to indicate when nodes
      are component of complexes
    * complexes are named as follows:
        * cp1&cp2&cp3
        * where cp1, cp2 and cp3 are the complex components
    * by the way, in this example of a 3 components complex, the following edges
      would be added as explain above:
        * cp1    membership    cp1&cp2&cp3
        * cp2    membership    cp1&cp2&cp3
        * cp3    membership    cp1&cp2&cp3

If -g/--geneSymbol is used:
    * kgml2sif attempts to name gene nodes according to a provided conversion
      table
    * the conversion table must be a 2 columns TSV file (see the file
      gene2symbol.tsv provided with kgml2sif)
    * unmatched KEGG IDs are left unchanged in the output SIF file and kgml2sif
      warns about them

Node name prefixes in KGML:
    * specify the node type
    * the name of nodes representing human genes is prefixed with "hsa:"
    * the name of nodes representing compounds is prefixed with "cpd:"
    * the name of nodes representing glycans (a subtype of compounds) is
      prefixed with "gl:"
    * the name of nodes representing drugs (a subtype of compounds) is prefixed
      with "dr:"

Full explanation of the KGML file format at https://www.kegg.jp/kegg/xml/docs/

For explanations about the SIF file format, see the readme file of kgml2sif.
""",
formatter_class=argparse.RawDescriptionHelpFormatter
)
parser.add_argument(
    "kgmlFile",
    type=str,
    metavar="<kgmlFile>",
    help="a KGML encoded KEGG human pathway",
    nargs="+"
)
parser.add_argument(
    "-g",
    "--geneSymbol",
    type=str,
    metavar="<tsvFile>",
    help="a conversion table for translating KEGG human gene IDs to gene symbols (see the file gene2symbol.tsv provided with kgml2sif)"
)
args=parser.parse_args()
#### GLOBALS ###################################################################
nodeTypes=[
    "gene",
    "compound",
    "group"
]
componentTypes=[
    "gene",
    "compound"
]
relationTypes=[
    "PPrel",
    "GErel",
    "PCrel",
    "ECrel"
]
gene2symbol={}
#### GENE2SYMBOL ###############################################################
if args.geneSymbol!=None:
    warnFile=list(os.path.splitext(args.geneSymbol))
    warnFile[1]=".warns.txt"
    warnFile="".join(warnFile)
    warns=[]
    for line in open(args.geneSymbol,"rt").read().splitlines():
        if line.count("\t")!=1:
            warns.append(line+": must be a 2 fields TSV")
        else:
            line=line.split("\t")
            ok=True
            for val in line:
                if val.strip()=="":
                    warns.append("\t".join(line)+": empty field")
                    ok=False
                    break
            if ok:
                if re.fullmatch(r"^hsa:[0-9]+$",line[0])==None:
                    warns.append("\t".join(line)+": invalid KEGG human gene ID")
                elif line[0] in gene2symbol.keys():
                    warns.append("\t".join(line)+": duplicated KEGG human gene ID")
                else:
                    gene2symbol[line[0]]=line[1]
    if len(warns)!=0:
        print("Warning: "+args.geneSymbol+": see "+warnFile)
        open(warnFile,"w").write("\n".join(warns)+"\n")
    if len(gene2symbol.keys())==0:
        print("Error: "+args.geneSymbol+": no valid lines")
        quit()
#### KGML FILES ################################################################
for kgmlFile in args.kgmlFile:
    #### OUT FILES #############################################################
    warnFile=list(os.path.splitext(kgmlFile))
    notFoundFile=list(os.path.splitext(kgmlFile))
    multiFile=list(os.path.splitext(kgmlFile))
    sifFile=list(os.path.splitext(kgmlFile))
    warnFile[1]=".warns.txt"
    notFoundFile[1]=".no_symbol.txt"
    multiFile[1]=".multiedges.sif"
    sifFile[1]=".sif"
    warnFile="".join(warnFile)
    notFoundFile="".join(notFoundFile)
    multiFile="".join(multiFile)
    sifFile="".join(sifFile)
    #### LOCALS ################################################################
    warns=[]
    notFound=[]
    nodes={}
    complexes={}
    edges={}
    memberships=[]
    regulars=[]
    multiEdges=[]
    sifNoDups=[]
    multiNoDups=[]
    root=xml.etree.ElementTree.parse(kgmlFile).getroot()
    #### NODES #################################################################
    for element in root.findall("entry"):
        type_=element.get("type")
        ID=element.get("id")
        name=element.get("name")
        if type_ in nodeTypes:
            if (ID==None) or (ID.strip()==""):
                warns.append(type_+": no ID")
            elif (ID in nodes.keys()) or (ID in complexes.keys()):
                warns.append(type_+": "+ID+": duplicated ID")
            elif (name==None) or (name.strip()==""):
                warns.append(type_+": "+ID+": no name")
            else:
                if type_=="group":
                    complexes[ID]=copy.deepcopy(element)
                else:
                    keggIDs=[]
                    ok=True
                    for keggID in name.split(" "):
                        if not any([
                            re.fullmatch(r"^hsa:[0-9]+$",keggID)!=None,
                            re.fullmatch(r"^cpd:C[0-9]+$",keggID)!=None,
                            re.fullmatch(r"^gl:G[0-9]+$",keggID)!=None,
                            re.fullmatch(r"^dr:D[0-9]+$",keggID)!=None
                        ]):
                            warns.append(type_+": "+ID+": "+keggID+": invalid name")
                            ok=False
                            break
                        elif keggID not in keggIDs:
                            keggIDs.append(keggID)
                    if ok:
                        if len(keggIDs)==0:
                            warns.append(type_+": "+ID+": "+name+": no valid name")
                        else:
                            nodes[ID]={
                                "type":type_,
                                "name":copy.deepcopy(keggIDs)
                            }
    #### GENE SYMBOLS ##########################################################
    if args.geneSymbol!=None:
        for ID in nodes.keys():
            if nodes[ID]["type"]=="gene":
                names=[]
                for i in range(len(nodes[ID]["name"])):
                    if nodes[ID]["name"][i] in gene2symbol.keys():
                        nodes[ID]["name"][i]=gene2symbol[nodes[ID]["name"][i]]
                    elif nodes[ID]["name"][i] not in notFound:
                        notFound.append(nodes[ID]["name"][i])
                for name in nodes[ID]["name"]:
                    if name not in names:
                        names.append(name)
                nodes[ID]["name"]=copy.deepcopy(names)
    #### COMPLEXES #############################################################
    for ID in complexes.keys():
        components=complexes[ID].findall("component")
        if len(components)==0:
            warns.append("group: "+ID+": no components")
        elif len(components)==1:
            warns.append("group: "+ID+": only one component")
        else:
            componentIDs=[]
            ok=True
            for component in components:
                componentID=component.get("id")
                if (componentID==None) or (componentID.strip()==""):
                    warns.append("group: "+ID+": component with no ID")
                    ok=False
                    break
                elif componentID not in nodes.keys():
                    warns.append("group: "+ID+": "+componentID+": unknown component")
                    ok=False
                    break
                elif nodes[componentID]["type"] not in componentTypes:
                    warns.append("group: "+ID+": "+componentID+": "+nodes[componentID]["type"]+": invalid component type")
                    ok=False
                    break
                elif componentID not in componentIDs:
                    componentIDs.append(componentID)
            if ok:
                if len(componentIDs)==0:
                    warns.append("group: "+ID+": no valid components")
                elif len(componentIDs)==1:
                    warns.append("group: "+ID+": only one valid component")
                else:
                    components=[]
                    cpxs=[]
                    nodes[ID]={
                        "type":"complex",
                        "name":[]
                    }
                    for componentID in componentIDs:
                        components.append(copy.deepcopy(nodes[componentID]["name"]))
                    m=len(components)
                    n=1
                    for component in components:
                        n*=len(component)
                    for i in range(n):
                        cpx=[]
                        for j in range(m):
                            cpx.append(None)
                        cpxs.append(copy.deepcopy(cpx))
                    for j in range(m):
                        q=len(components[j])
                        p=1
                        for i in range(j+1,m):
                            p*=len(components[i])
                        for i in range(n):
                            cpxs[i][j]=components[j][(i//p)%q]
                    for cpx in cpxs:
                        name="&".join(sorted(cpx))
                        if name not in nodes[ID]["name"]:
                            nodes[ID]["name"].append(name)
                        for node in cpx:
                            memberships.append([node,["membership"],name])
    #### RELATIONS #############################################################
    for relation in root.findall("relation"):
        type_=relation.get("type")
        entry1=relation.get("entry1")
        entry2=relation.get("entry2")
        subtypes=relation.findall("subtype")
        if type_ in relationTypes:
            if (entry1==None) or (entry1.strip()==""):
                warns.append(type_+": no source node")
            elif (entry2==None) or (entry2.strip()==""):
                warns.append(type_+": no target node")
            elif entry1 not in nodes.keys():
                warns.append(type_+": "+entry1+": unknown source node")
            elif entry2 not in nodes.keys():
                warns.append(type_+": "+entry2+": unknown target node")
            else:
                ok=True
                if len(subtypes)==0:
                    names=["unknown"]
                else:
                    names=[]
                    for subtype in subtypes:
                        name=subtype.get("name")
                        if (name==None) or (name.strip()==""):
                            warns.append(type_+": ("+entry1+","+entry2+"): empty name")
                            ok=False
                            break
                        elif name not in names:
                            names.append(name)
                if ok:
                    if len(names)==0:
                        warns.append(type_+": ("+entry1+","+entry2+"): no valid name")
                    else:
                        for source in nodes[entry1]["name"]:
                            for target in nodes[entry2]["name"]:
                                regulars.append([source,copy.deepcopy(names),target])
    #### MERGE REGULARS ########################################################
    for line in regulars:
        edges[line[0]]={}
    for line in regulars:
        edges[line[0]][line[2]]=[]
    for line in regulars:
        for relation in line[1]:
            if relation not in edges[line[0]][line[2]]:
                edges[line[0]][line[2]].append(relation)
    regulars=[]
    for source in edges.keys():
        for target in edges[source].keys():
            regulars.append([source,copy.deepcopy(edges[source][target]),target])
    #### MULTIEDGES ############################################################
    for line in memberships:
        if line[0] in edges.keys():
            if line[2] in edges[line[0]].keys():
                multiEdges.append(copy.deepcopy(line))
                multiEdges.append([line[0],copy.deepcopy(edges[line[0]][line[2]]),line[2]])
    #### SIF ###################################################################
    for i in range(len(memberships)):
        memberships[i][1]=",".join(sorted(memberships[i][1]))
        memberships[i]="\t".join(memberships[i])
    for i in range(len(regulars)):
        regulars[i][1]=",".join(sorted(regulars[i][1]))
        regulars[i]="\t".join(regulars[i])
    for i in range(len(multiEdges)):
        multiEdges[i][1]=",".join(sorted(multiEdges[i][1]))
        multiEdges[i]="\t".join(multiEdges[i])
    sif=copy.deepcopy(memberships)+copy.deepcopy(regulars)
    #### DUPLICATES ############################################################
    for line in sif:
        if line not in sifNoDups:
            sifNoDups.append(line)
    for line in multiEdges:
        if line not in multiNoDups:
            multiNoDups.append(line)
    sif=copy.deepcopy(sifNoDups)
    multiEdges=copy.deepcopy(multiNoDups)
    #### OUTPUTS ###############################################################
    if len(warns)!=0:
        print("Warning: "+kgmlFile+": see "+warnFile)
        open(warnFile,"w").write("\n".join(warns)+"\n")
    if len(notFound)!=0:
        print("Warning: "+sifFile+": contains unmatched KEGG IDs, see "+notFoundFile)
        open(notFoundFile,"w").write("\n".join(sorted(notFound))+"\n")
    if len(multiEdges)!=0:
        print("Warning: "+sifFile+": contains multiedges, see "+multiFile)
        open(multiFile,"w").write("\n".join(sorted(multiEdges))+"\n")
    if len(sif)==0:
        print("Error: "+kgmlFile+": no valid edges")
    else:
        open(sifFile,"w").write("\n".join(sorted(sif))+"\n")
