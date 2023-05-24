# Guide For Adding A New KEGG Pathway.
<ol>
<li> Log in to KEGG, select the pathway and download as in KGML format.
<li> Clone this Github repository: https://github.com/arnaudporet/kgml2sif. This is a KGML to simple interaction format tool.
<li> Run the tool with the following parameters: python kgml2sif.py -g conv/gene2symbol.tsv hsa00000.xml here we assume hsa00000.xml is the pathway name.
<li> The result is a table with intraction between proteins and other pathway compounds. You can select which ones you need for your pathway.
<li> In Pathway database each line represents either the input or the output of an interaction. The input or the output can a protein, a complex of a few proteins or a compound.
<li> Using Excel you need to provide the following table of interactions for this pathway:
	<ol>
		<li> column 1: pathway name.
		<li> column 2: pathway ID.
		<li> column 3: molecule type (i.e complex, protein, compound).
		<li> column 4: molecule names.
		<li> column 5: a unique complex or molecule number.
		<li> column 6: a comma-separated list of molecules names that are involved in this interaction. That is a bit different than the default KGML format. For example if we get:
			<ol>		
				KRAS	activation	ARAF
				KRAS	activation	BRAF
				KRAS	activation	PIK3CA
				KRAS	activation	PIK3CB
				KRAS	activation	PIK3CD
				KRAS	activation	PIK3R1
				KRAS	activation	PIK3R2
				
				We need to group the molecules to one list like the following:
				- path_name \t path_ID \t mol_type \t KRAS \t KRAS \t mol_num \t ARAF, BRAF, PIK3CA, PIK3CB, PIK3CD, PIK3R1, PIK3R2
			</ol>

		<li> column 7: optional: you can add 'active' in case of an active molecule. This will add a '+' sign in the KGML export parser to signal an active molecule.
		<li> column 8: add the source database, in this case, 'KEGG'.
		<li> column 9: ignored.
		<li> column 10: molecule role. Can be: input, output or inhibitor (only in case it negatively affects the interaction).
		<li> column 11: unique interaction ID.
		<li> column 12: interaction type: for example modification, degradation, translocation. Used only in graphics.
	</ol>
<li> After the table is prepared it needs to be concatenated to the pathologist.db.txt file in the data folder.
</ol>