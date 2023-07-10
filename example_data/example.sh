
# dot plot
# A
soi dotplot -s Populus_trichocarpa-Salix_dunnii.collinearity.gz \
	-g Populus_trichocarpa-Salix_dunnii.gff.gz -c Populus_trichocarpa-Salix_dunnii.ctl  \
	--kaks Populus_trichocarpa-Salix_dunnii.collinearity.ks.gz \
	--xlabel '$Populus\ trichocarpa$' --ylabel '$Salix\ dunnii$' \
	--ks-hist --max-ks 1.5 -o Populus_trichocarpa-Salix_dunnii     \
	--plot-ploidy --gene-axis --number-plots  
# B
soi dotplot -s Populus_trichocarpa-Salix_dunnii.orthologs.gz    \
	-g Populus_trichocarpa-Salix_dunnii.gff.gz -c Populus_trichocarpa-Salix_dunnii.ctl  \
	--kaks Populus_trichocarpa-Salix_dunnii.collinearity.ks.gz \
	--xlabel '$Populus\ trichocarpa$' --ylabel '$Salix\ dunnii$'  \
	--ks-hist --max-ks 1.5 -o Populus_trichocarpa-Salix_dunnii.o   \
	--plot-ploidy --gene-axis --number-plots  \
	--homology	# homology input
# C
soi dotplot -s Populus_trichocarpa-Salix_dunnii.collinearity.gz \
	-g Populus_trichocarpa-Salix_dunnii.gff.gz -c Populus_trichocarpa-Salix_dunnii.ctl  \
	--xlabel '$Populus\ trichocarpa$' --ylabel '$Salix\ dunnii$' \
	--ks-hist --max-ks 1 -o Populus_trichocarpa-Salix_dunnii.io    \
	--plot-ploidy --gene-axis --number-plots \
	--ofdir OrthoFinder/OrthoFinder/Results_*/ --of-color	# coloring by Orthology Index
# D
soi dotplot -s Populus_trichocarpa-Salix_dunnii.collinearity.gz \
	-g Populus_trichocarpa-Salix_dunnii.gff.gz -c Populus_trichocarpa-Salix_dunnii.ctl  \
	--kaks Populus_trichocarpa-Salix_dunnii.collinearity.ks.gz \
	--xlabel '$Populus\ trichocarpa$' --ylabel '$Salix\ dunnii$' \
	--ks-hist --max-ks 1.5 -o Populus_trichocarpa-Salix_dunnii.io  \
	--plot-ploidy --gene-axis --number-plots \
	--ofdir OrthoFinder/OrthoFinder/Results_*/ --of-ratio 0.6	# filtering by Orthology Index


# filter
soi filter -s Populus_trichocarpa-Salix_dunnii.collinearity.gz -o OrthoFinder/OrthoFinder/Results_*/ \
	-c 0.6 > Populus_trichocarpa-Salix_dunnii.collinearity.ortho

