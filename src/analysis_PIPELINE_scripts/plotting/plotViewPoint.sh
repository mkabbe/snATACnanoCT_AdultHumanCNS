
#hicPlotViewpoint --matrix P27454_hOPC_REdeep/P27454_REdeep_hOPC.mapped.pairs.contact-map.10000.cool \
#	--region "chr7:26460996-27860500" \
#		-o hoxA_tad_viewpoint.pdf \
#			-rp "chr7:27156350-27161919" \
#				--dpi 600


#hicPlotViewpoint --matrix P27454_hOPC_REdeep/P27454_REdeep_hOPC.mapped.pairs.contact-map.10000.cool \
#	--region "chr2:176000000-177200000" \
#		-o LINC01116_viewpoint.pdf \
#			-i LINC01116_viewpoint.bedgraph \
#			-rp "chr2:176634964-176639984" \
#				--dpi 600

#hicPlotViewpoint --matrix cools/hOPC_10kb.cool \
#hicPlotViewpoint --matrix P27454_hOPC_REdeep/P27454_REdeep_hOPC.mapped.pairs.contact-map.10000.cool \
#	--region "chr2:175600000-176700000" \
#		-o MIR10B_viewpoint.pdf \
#			-i MIR10B_viewpoint.bedgraph \
#			-rp "chr2:176150303-176150412" \
#				--dpi 600

hicPlotViewpoint --matrix cools/hOPC_5kb.cool \
        --region "chr7:26200000-28500000" \
                -o plots/HOXA_viewpoint_hOPC.pdf \
                        -i HOXA_viewpoint_hOPC.bedgraph \
                                -rp "chr7:27091875-27207580" \
                                        --dpi 600

hicPlotViewpoint --matrix cools/bCell3_5kb.cool \
	--region "chr7:26200000-28500000" \
		-o plots/HOXA_viewpoint_bCell.pdf \
			-i HOXA_viewpoint_bCell.bedgraph \
				-rp "chr7:27091875-27207580" \
					--dpi 600
