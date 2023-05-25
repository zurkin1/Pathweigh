# Guide For Adding A New KEGG Pathway.
<ol>
<li> Log in to KEGG, select the pathway and download as in KGML format.
<li> Clone this Github repository: https://github.com/arnaudporet/kgml2sif. This is a KGML to simple interaction format tool.
<li> Run the tool with the following parameters: python kgml2sif.py -g conv/gene2symbol.tsv hsa00000.xml here we assume hsa00000.xml is the pathway name.
<li> The result is a table with intraction between proteins and other pathway compounds. You can select which ones you need for your pathway.
<li> In Pathway database each line represents a part of an interaction, either an input or an output. The input or the output can be a protein, a complex of a few proteins or a compound.
<li> Using Excel you need to provide the following table of interactions for this pathway:
	<ol>
		<li> column 1: pathway name.
		<li> column 2: pathway ID.
		<li> column 3: molecule type (i.e protein, compound). compound will always get probabilty 1 in calculations.
		<li> column 4: molecule name. Used only in graphics. Must be repeated in column 6.
		<li> column 5: a unique molecule number.
		<li> column 6: a comma-separated list of molecules names that are involved in this interaction. That is a bit different than the default KGML format. For example if we have:
			<ol>
				<li>	KRAS	activation	ARAF
				<li>	KRAS	activation	BRAF
				<li>	KRAS	activation	PIK3CA
				<li>	KRAS	activation	PIK3CB
				<li>	KRAS	activation	PIK3CD
				<li>	KRAS	activation	PIK3R1
				<li>	KRAS	activation	PIK3R2	 We need to group the molecules to one list like the following. So two lines describe an interaction.
				<li>    path_name \t path_ID \t mol_type \t KRAS \t mol_num \t 	 
				<li>    path_name \t path_ID \t mol_type \t ARAF \t mol_num \t BRAF, PIK3CA, PIK3CB, PIK3CD, PIK3R1, PIK3R2
			</ol>
		<li> column 7: optional: you can add 'active' in case of an active molecule. This will add a '+' sign in the KGML export parser to signal an active molecule.
		<li> column 8: add the source database, in this case, 'KEGG'.
		<li> column 9: ignored.
		<li> column 10: The molecule role in column 4. Can be: input, output or inhibitor (only in case it negatively affects the interaction).
		<li> column 11: unique interaction ID.
		<li> column 12: interaction type. Can be any detailed string, for example modification, degradation, translocation. Used only in graphics.
	</ol>
<li> After the table is prepared it needs to be renamed to pathologist.db.txt file in the data folder. Don't concatenate to the old database file, use them separately.
<li> An example of a file in this format is available in

###[pathologist.db.txt](keggpathologist.db.txt)
</ol>