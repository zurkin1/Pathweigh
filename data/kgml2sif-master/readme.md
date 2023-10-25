# Converting KGML encoded KEGG human pathways to the SIF file format

Copyright (C) 2019-2020 [Arnaud Poret](https://github.com/arnaudporet)

This work is licensed under the [GNU General Public License](https://www.gnu.org/licenses/gpl.html).

To view a copy of this license, visit https://www.gnu.org/licenses/gpl.html.

## kgml2sif

### Requirements

* [Python 3](https://www.python.org)
* a unix based operating system (ex: [Linux](https://en.wikipedia.org/wiki/Linux), [macOS](https://en.wikipedia.org/wiki/MacOS))

### Usage

Ensure that `kgml2sif` is executable:

```sh
chmod ugo+rx kgml2sif
```

Usage: `kgml2sif [-h] [-g <tsvFile>] <kgmlFile> [<kgmlFile> ...]`

Positional arguments:

* `<kgmlFile>`: a KGML encoded KEGG human pathway

Optional arguments:

* `-g <tsvFile>`, `--geneSymbol <tsvFile>`: a conversion table for translating KEGG human gene IDs to gene symbols (see the file `gene2symbol.tsv` provided with kgml2sif in the `conv` folder)
* `-h`, `--help`: show help

### Help

Convert KGML encoded KEGG human pathways to the SIF file format.

Currently:

* the following node types are considered:
    * `gene` (in KGML, genes also stand for the corresponding proteins/gene products)
    * `compound`
    * `group` (in KGML, groups stand for complexes)
* the following relation types are considered (in KGML, relations stand for edges):
    * `PPrel`: protein-protein relations
    * `GErel`: protein-gene relations (i.e. gene expression regulations)
    * `PCrel`: protein-compound relations
    * `ECrel`: enzyme-enzyme relations sharing a common compound (as product for the first enzyme and as substrate for the second one)

In the output SIF file:

* edges for which the name is missing are automatically named `unknown`
* multiedges, if any, are merged if possible, otherwise they are left inside the output SIF file and kgml2sif warns about them
* the non KGML `membership` edges are added in order to indicate when nodes are component of complexes
* complexes are named as follows:
    * `cp1&cp2&cp3`
    * where `cp1`, `cp2` and `cp3` are the complex components
* by the way, in this example of a 3 components complex, the following edges would be added as explain above:
    * `cp1    membership    cp1&cp2&cp3`
    * `cp2    membership    cp1&cp2&cp3`
    * `cp3    membership    cp1&cp2&cp3`

If `-g/--geneSymbol` is used:

* kgml2sif attempts to name gene nodes according to a provided conversion table
* the conversion table must be a 2 columns TSV file (see the file `gene2symbol.tsv` provided with kgml2sif in the `conv` folder)
* unmatched KEGG IDs are left unchanged in the output SIF file and kgml2sif warns about them

Node name prefixes in KGML:

* specify the node type
* the name of nodes representing human genes is prefixed with `hsa:`
* the name of nodes representing compounds is prefixed with `cpd:`
* the name of nodes representing glycans (a subtype of compounds) is prefixed with `gl:`
* the name of nodes representing drugs (a subtype of compounds) is prefixed with `dr:`

Full explanation of the KGML file format at https://www.kegg.jp/kegg/xml/docs/.

For explanations about the SIF file format, see at the end of this file.

## Examples

These examples come from downloaded human [KEGG pathways](https://www.genome.jp/kegg/pathway.html).

The SIF files produced by kgml2sif from these examples are also in the SVG file format for visualization purpose.

### Apoptosis

```sh
./kgml2sif -g conv/gene2symbol.tsv examples/Apoptosis/Apoptosis.xml
```

### ErbB signaling pathway

```sh
./kgml2sif -g conv/gene2symbol.tsv examples/ErbB_signaling_pathway/ErbB_signaling_pathway.xml
```

### Insulin signaling pathway

```sh
./kgml2sif -g conv/gene2symbol.tsv examples/Insulin_signaling_pathway/Insulin_signaling_pathway.xml
```

### p53 signaling pathway

```sh
./kgml2sif -g conv/gene2symbol.tsv examples/p53_signaling_pathway/p53_signaling_pathway.xml
```

### TNF signaling pathway

```sh
./kgml2sif -g conv/gene2symbol.tsv examples/TNF_signaling_pathway/TNF_signaling_pathway.xml
```

## The KGML file format

KGML stands for KEGG Markup Language. It is an [XML](https://www.w3.org/XML/) representation of [KEGG pathways](https://www.genome.jp/kegg/pathway.html) and is a file format in which KEGG pathways can be downloaded, either from the KEGG web site or using the [KEGG API](https://www.kegg.jp/kegg/rest/keggapi.html).

For a complete explanation of the KGML file format, see https://www.kegg.jp/kegg/xml/docs/.

## The SIF file format

In a SIF file encoding a network, each line encodes an edge as follows:

```
source \t interaction \t target
```

Note that the field separator is the tabulation: the SIF file format is the tab separated values format (TSV) with exactly 3 columns.

For example, the edge representing the activation of RAF1 by HRAS is a line of a SIF file encoded as follows:

```
HRAS \t activation \t RAF1
```
