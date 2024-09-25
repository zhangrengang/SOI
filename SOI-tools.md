
#### Macro-synteny phylogeny ####
Before run the pipeline:
1. install the lasest verion of [SOI](https://github.com/zhangrengang/orthoindex#installation) by `python3 setup.py install`, 
and have a check: `which soi-syn`;
2. completed the [example pipeline](https://github.com/zhangrengang/evolution_example) to get the orthologous synteny file `collinearity.ortho`;
3. prepare file `species.config` to set the expected subgenome numbers for the targeted species (TAB seperated):
```
Vitis_vinifera  1
Aralia_elata    2
Centella_asiatica       2
```
Then run the pipeline:
```
cd phylogenomics
soi-syn anchor_trees collinearity.ortho species.config ../pep.faa ../all_species_gene.gff output_dir
```
After the pipeline completed, you can find tree files in the `output_dir`:
```
OG*treefile		# gene tree file for each anchor gene
chr*treefile	# macro-synteny tree file by concatenating the anchor genes from the same chromosome set
CHR*treefile	# macro-synteny tree file, but allowing gene missing
```
For example, `CHR_Ae11-Ae15-Ca2-Ca7-Vv1_143_283.concat.treefile`, `Ae11-Ae15-Ca2-Ca7-Vv1` is the chromosome set,
`143` is the number of anchor genes, and `283` is the number of all syntenic genes allowing missing.

