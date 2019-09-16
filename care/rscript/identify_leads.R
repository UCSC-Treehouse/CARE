

# usage example (after running function definition below ) 
sampleID="TH06_0646_S01"; sampleDataDir="~/Documents/Dropbox/ucsc/operations/tumor board prep/POG-BC June 8/TH06_0646_S01_scratch"; pathwayAnnotationFileName = "../review of outliers/curatedPathwaysContainingFDA_druggableGenes.txt"; FDA_DruggableGenesByCategoryFile="../review of outliers/genesByCategory.txt"

getDruggableOutliersInEnrichedPathways(sampleID, sampleDataDir, pathwayAnnotationFileName, FDA_DruggableGenesByCategoryFile)



getDruggableOutliersInEnrichedPathways<-function(sampleID, sampleDataDir, pathwayAnnotationFileName, FDA_DruggableGenesByCategoryFile){

	library(tidyverse)
	
	setwd(sampleDataDir)
	pathwayAnnotation=read_tsv(pathwayAnnotationFileName, col_names=FALSE)
	colnames(pathwayAnnotation)=c("geneSet", "anno")
	
	FDA_DruggableGenesByCategory= read_tsv(FDA_DruggableGenesByCategoryFile)



	# replace previous script
	multiSampleData=read_tsv(multiSampleDataFile)
	
	###
	### FDA druggable genes
	###
	
	FDA_DruggableGenesByCategory= read_tsv(FDA_DruggableGenesByCategoryFile)
	
	FDA_DruggableGenesByCategory$group=factor(FDA_DruggableGenesByCategory$group, levels=unique(FDA_DruggableGenesByCategory$group))
	
	s= sampleID
	
	#setwd(file.path("~/Documents/Dropbox/ucsc/operations/tumor board prep/POG-BC June 8/", s))
	
	enrichedPathways =read_tsv(paste0(s, " GeneSetAggregation.txt"))
	enrichedPathways $enrichedOutsideTopFivePercent=rowSums(enrichedPathways[,2:4])>0

#	enrichedPathways$tumorType= multiSampleData $`Primary group` [match(s, multiSampleData $`Sample THID`)]
#	enrichedPathways$tumorType= importedDataList[["importedSamplesDF"]]$`Primary group` [match(s, importedDataList[["importedSamplesDF"]]$`Sample THID`)]
	pathwayDetailsRaw=read_tsv(paste0(s, " GeneSetDetailsPerList.txt"))
	pathwayDetailsRaw $comm_up_kkText=paste0(pathwayDetailsRaw $`comm_up_N_Genes_in_Overlap_(k)`, "/", pathwayDetailsRaw $`N_Genes_in_Gene_Set_(K)`)
	pathwayDetailsRaw $pc_up_kkText=paste0(pathwayDetailsRaw $`pc_up_N_Genes_in_Overlap_(k)`, "/", pathwayDetailsRaw $`N_Genes_in_Gene_Set_(K)`)
	pathwayDetailsRaw $pd_up_kkText=paste0(pathwayDetailsRaw $`pd_up_N_Genes_in_Overlap_(k)`, "/", pathwayDetailsRaw $`N_Genes_in_Gene_Set_(K)`)
	
	enrichedPathways$comm_up_kkText= pathwayDetailsRaw $comm_up_kkText[match(enrichedPathways$GeneSet, pathwayDetailsRaw$GeneSetName)]
	enrichedPathways$pc_up_kkText= pathwayDetailsRaw $pc_up_kkText[match(enrichedPathways$GeneSet, pathwayDetailsRaw$GeneSetName)]
	enrichedPathways$pd_up_kkText= pathwayDetailsRaw $pd_up_kkText[match(enrichedPathways$GeneSet, pathwayDetailsRaw$GeneSetName)]
	
	
	kkcols=which(grepl("kkText", colnames(enrichedPathways)))
	enrichedPathways $textOfEnrichmentList<-NA
	# generate textOfEnrichmentList
	for (i in 1:nrow(enrichedPathways)){
		# i=1
		whichSets=which(enrichedPathways[i,2:4]==1)
		t2=paste0( colnames(enrichedPathways)[2:4][whichSets], " (", enrichedPathways[i, kkcols][whichSets], ")")
		enrichedPathways $textOfEnrichmentList [i]=gsub(", enriched_in_", " and ", paste(t2, collapse=", "))
	}
	# t1= enrichedPathways[,2:4]==1
	# t2=apply(t1, 1, function(x) paste( colnames(enrichedPathways)[2:4][x] )
	# t3=lapply(t2, function(x){
		# gsub(", enriched_in_", " and ", paste(x, collapse=", "))
		# })
	# enrichedPathways $textOfEnrichmentList=unlist(t3)
	
	enrichedPathways $geneSetAndTextOfEnrichmentList=paste0(enrichedPathways $GeneSet , " (", enrichedPathways $textOfEnrichmentList, ")")
	
	enrichedPathways$FDA_druggableGenesInGeneSet= lapply(enrichedPathways$allMemberGenesInThLists, function(x) strsplit(x, ", ")[[1]])
	
	enrichedPathways$countOfFDA_druggableGenesInGeneSet =unlist(lapply(enrichedPathways $FDA_druggableGenesInGeneSet, length))
	
	pathwayAnnotation$anno[is.na(pathwayAnnotation$anno)]="geneSets_priority1"
	
	enrichedPathways$anno= pathwayAnnotation$anno[match(enrichedPathways$GeneSet, pathwayAnnotation$geneSet)]
		
	pc_upBySample=read_tsv(paste0("genes_", s, "_pc_up"))
	pd_upBySample=read_tsv(paste0("genes_", s, "_pd_up"))
	
	colnames(pd_upBySample)[3]=colnames(pc_upBySample)[3]="median"
	
	upOutliers=data.frame(rbind(cbind(pc_upBySample, outlierComparisonGroup="pc"), cbind(pd_upBySample, outlierComparisonGroup="pd")))
	
	colnames(upOutliers)=c("gene", "log2TPMp1", "groupMedian", "outlierComparisonGroup")
	
	upOutliers$category= FDA_DruggableGenesByCategory$group[match(upOutliers$gene, FDA_DruggableGenesByCategory$gene)]
	upOutliers$hasFdaApprovedTargetedDrug=!is.na(upOutliers$category)
	
	expectedColCount=length(unique(FDA_DruggableGenesByCategory$group ))+1 # one column for sample name
	
	theseOutliers=subset(upOutliers,  gene %in% FDA_DruggableGenesByCategory$gene )
	theseEnrichedPathways=subset(enrichedPathways,  enrichedOutsideTopFivePercent & countOfFDA_druggableGenesInGeneSet>0)
	
	###############
	theseOutliers=subset(upOutliers,  gene %in% FDA_DruggableGenesByCategory$gene )
	theseEnrichedPathways=subset(enrichedPathways,  enrichedOutsideTopFivePercent & countOfFDA_druggableGenesInGeneSet>0)
	
	#
	# assemble FDA_DruggableGenesPerEnrichedDruggablePathway
	#

	FDA_DruggableGenesPerEnrichedDruggablePathway<-NULL
	outputDF<-NULL
	
	if (nrow(theseEnrichedPathways)>0 & nrow(theseOutliers)>0 ) {		
		theseOutliersByList=aggregate(outlierComparisonGroup ~ gene, theseOutliers, paste, collapse=" and ") # run above
		for (j in 1:nrow(theseEnrichedPathways)){	
			p=theseEnrichedPathways $GeneSet[j]
			thisPathwayInfo=subset(theseEnrichedPathways, GeneSet==p)
			genesInPathway=		unlist(thisPathwayInfo $FDA_druggableGenesInGeneSet[unlist(lapply(thisPathwayInfo $FDA_druggableGenesInGeneSet, length))>0])
			theseFDADruggableGenesInPathway=subset(theseOutliersByList, gene %in% genesInPathway)
			# add FDA gene set, e.g. PI3K/AKT/mTOR from theseOutliers
	
			thisPathwayDescription=gsub("enriched_in_", "", paste0(theseEnrichedPathways $GeneSet[j], " (", theseEnrichedPathways $textOfEnrichmentList[j], ")"))
	
			if (nrow(theseFDADruggableGenesInPathway)>0){
				theseFDADruggableGenesInPathway$geneDesc=paste0(theseFDADruggableGenesInPathway $gene, " (", theseFDADruggableGenesInPathway$outlierComparisonGroup, ")")
				FDA_DruggableGenesPerEnrichedDruggablePathway <-rbind(FDA_DruggableGenesPerEnrichedDruggablePathway,  data.frame(THid= sampleID, assay= thisPathwayInfo$anno, results= thisPathwayDescription, details=paste("druggable:", paste0(theseFDADruggableGenesInPathway$geneDesc, collapse=", ")), details2=paste("overlap_genes_in_pc_up", paste(thisPathwayInfo$overlap_genes_in_pc_up, "and pd_up", thisPathwayInfo$overlap_genes_in_pd_up ))))
			}
		}
	} else {
			FDA_DruggableGenesPerEnrichedDruggablePathway <-rbind(FDA_DruggableGenesPerEnrichedDruggablePathway,  data.frame(THid= sampleID, assay="geneSets", results= "no enriched druggable geneSets", details="", details2=""))
	}
	
	##########################################
	### Note FDA Druggable up outliers
	###################################
	
	if (nrow(theseOutliers)>0){
		allDruggableUpOutliers <-data.frame(THid= sampleID, assay="druggableUpOutlier", theseOutliers %>% group_by(gene) %>% summarize (outlierGroup=paste(outlierComparisonGroup, collapse=", ")))
		colnames(allDruggableUpOutliers)[3:4]=c("results", "details")
	} else {
		allDruggableUpOutliers <-data.frame(THid= sampleID, assay="druggableUpOutlier", results ="No druggableUpOutliers", details="")
	}
	allDruggableUpOutliers$details2=""
	
	
	outputDF=rbind(FDA_DruggableGenesPerEnrichedDruggablePathway, allDruggableUpOutliers)
	assayOrder=c("druggableUpOutlier", "geneSets_priority1", "geneSets_broadCancer", "geneSets_nonCancer" )
	outputDF $assay =factor(outputDF $assay, levels=unique(c(assayOrder, outputDF $assay)))
	outputDF <-arrange(outputDF, assay)
	
	
	write_tsv(outputDF, paste0("automatedLeadsIdentified.tsv"))
	
	
}