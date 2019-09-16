
# version notes
# v10 - fixed column names; now you can distinguish pc_up from pd_up
# v11 - add overview tab to excel output <- expected to change 
# v12 - fixed duplication of data in two columns: pc-dgidb and pcd-dgidb 
# v13 - added "unique" to fix duplication of data when retreieving drug interactions, e.g. for pancanUp_DGIdb
# v14 - add ranking to drugs
# v15 - export ranking to drugs in text as well as xl
# v16 - adds overlap_genes_in_comm_up	overlap_genes_in_pc_up	overlap_genes_in_pd_up	overlap_genes_in_top5 to GeneSetAggregation

# Libraries 
library(reshape2)
library(data.table)
# install.packages("openxlsx", dependencies=TRUE) # Use this if package is not installed
# library(openxlsx)
# Sys.setenv(R_ZIPCMD= "/usr/bin/zip")

# Function - combine gsea & dgidb into an excel worksheet
# Dir must not contain gsea_ and dgidb_ files that don't pertain to the sample

combine.gsea.dgidb<-function(sampleID, dataDir, pathwayFile){
    
    print(paste("Processing sample: ", sampleID)) # ETK DEBUG

    
    gseaEmptyFileSizeMaximum=175
    dgiDbEmptyFileSize=1

    ###
    ### setup
    ###

    options(stringsAsFactors=FALSE)

    read.txt<-function (file="", sep="\t", header=TRUE,row.names=NULL, ...){
        print(paste("reading file ", file)) # ETK DEBUG
        read.table(file, sep= sep, header= header, row.names= row.names, ...)
        }

    write.txt<-function (x, file="", quote=FALSE, sep="\t", row.names=FALSE, ...){
        write.table(x, file, quote=quote, sep=sep, row.names=row.names, ...)
        }

    ###
    ### load generic pathways
    ###
    #allPathwaysFileRaw="msigdb.v5.2.symbols.gmt"
    allPathwaysFileRaw=pathwayFile

    #setwd(pathwayFileLoc)

    # Don't  create processed allPathwaysFile  
    # allPathwaysFileProcessed=gsub(".gmt$", ".matrix.txt.gz", allPathwaysFileRaw)
    # if (file.exists(allPathwaysFileProcessed)){
    #     allPathwaysDF =fread(paste("gzip -d --stdout", allPathwaysFileProcessed), data.table = FALSE)
    # } else {
    allPathwaysRaw=scan(allPathwaysFileRaw, what="list", sep="\n")
    maxLen=0
    allPathwaysMatrixList=lapply(allPathwaysRaw, function(x) { y1=strsplit(x, "\t")[[1]]; y2= y1[3:length(y1)];matrix(data=c(rep(y1[1], length=length(y2)), y2), ncol=2, byrow=FALSE)})
    allPathwaysDF <-data.frame(do.call("rbind", allPathwaysMatrixList))
    colnames(allPathwaysDF)=c("GeneSet", "Gene")
    #    allPathwaysConn <- gzfile(allPathwaysFileProcessed)
    #     write.txt(allPathwaysDF, allPathwaysConn)
    # }	
    setwd(dataDir)


    ###
    ### load genes in lists, e.g. pancan, panDisease and top 5pct
    ###
    geneMeta=data.frame(fn=list.files(,paste0("genes_", sampleID, "*")))
    # currently, this imports the following files: genes_[sampleID]_comm_up, genes_[sampleID]_pc_down, genes_[sampleID]_pc_up, genes_[sampleID]_pd_down, genes_[sampleID]_pd_up, genes_[sampleID]_top5
    geneMeta $sampleName= sampleID# lgCommonSubstring(geneMeta $fileID)
    geneMeta $dataTag=gsub(paste0("^.*",sampleID, "_"), "", geneMeta $fn)
    geneMeta=subset(geneMeta, ! grepl(".xlsx$", fn))

    # make matrix of genes v lists
    geneDataRaw=lapply(geneMeta $fn, read.txt) #("genes_SRR1988322b")
    geneNameAndTag=Reduce("rbind",lapply(1:nrow(geneMeta), function(x) data.frame(gene= geneDataRaw[[x]][,1], sampleName= geneMeta $sampleName[x] ,dataTag = geneMeta $dataTag[x])))
    geneByThGeneList=dcast(geneNameAndTag, gene ~ dataTag, value.var="sampleName", fun.aggregate=length)
    oldnames=colnames(geneByThGeneList)[2:ncol(geneByThGeneList)]
    colnames(geneByThGeneList)=c("gene", paste0("inSet_", oldnames))
    geneByThGeneList$thListsContainingGene=rowSums(geneByThGeneList[2:ncol(geneByThGeneList)])


    ###
    ### identify which genes are druggable
    ###

    # get metadata from dgidb files
    dgMeta=data.frame(fn=list.files(,"dgidb_."))
    dgMeta$sampleName= sampleID
    dgMeta$dataTag=gsub(paste0("^.*",sampleID, "_"), "", dgMeta $fn)
    dgMeta$fileInfo=file.info(dgMeta $fn)
    dgMeta =subset(dgMeta, ! grepl(".xlsx$", fn))

    dgMeta $empty = dgMeta $fileInfo$size== dgiDbEmptyFileSize
    dgMetaNonEmpty=subset(dgMeta, !empty)


    # make matrix of genes v lists
    dgDataRaw=lapply(dgMetaNonEmpty $fn, read.txt, header=F) #("dgidb_export_2017-01-03.SRR1988322c.tsv")

    dgDataInteractions=data.frame(rawInteraction=unlist(lapply(dgDataRaw, function(x)  x$V1)))
    dgDataInteractions$gene=gsub(" .* .*", "", dgDataInteractions $rawInteraction)
    dgDataInteractions$drug=gsub("^.*and ", "", dgDataInteractions $rawInteraction)
    dgDataGeneNames=lapply(dgDataRaw, function(x) unique(gsub(" .* .*", "", x$V1)))
    dgDataGeneNameAndTag=Reduce("rbind",lapply(1:nrow(dgMetaNonEmpty), function(x) data.frame(druggableGene=dgDataGeneNames[[x]], sampleName=dgMeta$sampleName[x] ,dataTag = dgMetaNonEmpty $dataTag[x])))
    druggableGeneByThGeneList=dcast(dgDataGeneNameAndTag, druggableGene ~ dataTag, value.var="sampleName", fun.aggregate=length)
    druggableGeneByThGeneList $thListsContainingGene =rowSums(druggableGeneByThGeneList[2:ncol(druggableGeneByThGeneList)])
    druggableGeneByThGeneList = druggableGeneByThGeneList[order(druggableGeneByThGeneList $thListsContainingGene, decreasing=TRUE),]
    # write.txt(druggableGeneByThGeneList, file=paste(dgMeta$sampleName[1], "druggableGeneAggregation.txt")) ## moved to after

    geneByThGeneList$druggableGene= geneByThGeneList$gene %in% druggableGeneByThGeneList$druggableGene

    geneByThGeneList2= geneByThGeneList

    colsToCheck=2:7
    setsToConsider=colnames(geneByThGeneList2 ) [colsToCheck]

    geneByThGeneList2$setCombo=apply(geneByThGeneList2[, colsToCheck], 1, function(x) gsub("inSet_", "", paste(setsToConsider[as.logical(x)], collapse=",")))
    # x= geneByThGeneList2[21,]

    geneByThGeneList2$setCombo[geneByThGeneList2$druggableGene]=paste0(geneByThGeneList2$setCombo[geneByThGeneList2$druggableGene], ",druggable")

    write.txt(geneByThGeneList, paste(dgMeta$sampleName[1], "allGeneAggregation.txt"))


    ###
    ### analyze gene sets enriched pathways
    ###

    # get metadata from gsea files
    gseaMeta=data.frame(fn=list.files(,"gsea_.*"))
    gseaMeta$sampleName=sampleID # lgCommonSubstring(gseaMeta$fileID)
    gseaMeta$dataTag=gsub(paste0("^.*",sampleID, "_"), "", gseaMeta $fn)
    gseaMeta =subset(gseaMeta, ! grepl(".xlsx$", fn))


    gseaMeta $fileInfo=file.info(gseaMeta $fn)
    gseaMeta $empty = gseaMeta $fileInfo$size <= gseaEmptyFileSizeMaximum
    gseaMetaNonEmpty=subset(gseaMeta, !empty)



    gseaGeneSetListRaw <-list()

    ### identify pathways enriched per Th Gene List
    for (i in 1:nrow(gseaMetaNonEmpty)){
        # i=1
        #
        # pull locations out of multi-table gsea file
        #
        allGseaInfoRaw=scan(gseaMetaNonEmpty $fn[i], what="list", sep="\n", blank.lines.skip=FALSE)
        firstLineOfGeneSetTable=grep("^Gene Set Name", allGseaInfoRaw)
        afterEndOfGeneSetTable=grep("Gene/Gene Set Overlap Matrix", allGseaInfoRaw)
        allBlankLines=which(allGseaInfoRaw =="")
        lastLineOfGeneSetTable =sort(allBlankLines[allBlankLines<afterEndOfGeneSetTable],decreasing=TRUE)[2]-1
        firstLineOfGeneSet_GeneMatrix=grep("^Entrez Gene Id", allGseaInfoRaw)

        #
        # pull gene set name table out of multi-table gsea file
        #
        geneSetTable=read.txt(gseaMetaNonEmpty $fn[i], fill=T, comment.char="", quote="", skip= firstLineOfGeneSetTable-1, nrows= lastLineOfGeneSetTable -firstLineOfGeneSetTable)
        # geneSetTable= allGseaInfo[(1+firstLineOfGeneSetTable): lastLineOfGeneSetTable,1:7]
        colnames(geneSetTable)=gsub(" ", "_", c("GeneSetName", "N Genes in Gene Set (K)", "Description", "N Genes in Overlap (k)", "k/K", "p-value", "FDR q-value")) # dput(as.character(allGseaInfo[firstLineOfGeneSetTable,1:7]))

        gseaGeneSetListRaw[[i]]= geneSetTable
        names(gseaGeneSetListRaw[i])= gseaMetaNonEmpty $fileID[i]
    }	
	
	
	###
	### create a table listing all reported gene sets and identify which ThGeneLists they're enriched in
	###
	
	names(gseaGeneSetListRaw)= gseaMetaNonEmpty$dataTag
	
	gseaEnrichedGeneSetsList=lapply(gseaGeneSetListRaw, function(x) unique(x$GeneSetName))
	
	gseaEnrichedGeneSetsByThGeneList=Reduce("rbind",lapply(1:nrow(gseaMetaNonEmpty), function(x) data.frame(GeneSet= gseaEnrichedGeneSetsList[[x]], sampleName= gseaMetaNonEmpty $sampleName[x] ,dataTag = gseaMetaNonEmpty $dataTag[x])))
	
	gseaEnrichedGeneSetsByThGeneList$dataTag=paste0("enriched_in_", gseaEnrichedGeneSetsByThGeneList$dataTag)
	
	enrichedGeneSet=dcast(gseaEnrichedGeneSetsByThGeneList[,c("GeneSet", "sampleName", "dataTag")], GeneSet ~ dataTag, value.var="sampleName", fun.aggregate=length)
	
	# add cols when nothing in list is enriched, e.g. "enriched_in_pc_up"
	expectedThGeneSetsWithGeneSetAnalysis=c("comm_up", "pc_up", "pd_up", "top5")
	enrichedColsToAdd=! paste0("enriched_in_", expectedThGeneSetsWithGeneSetAnalysis)  %in% colnames(enrichedGeneSet)
	if (sum(enrichedColsToAdd)>0){
		emptyEnrichedToAdd= data.frame(matrix(, ncol=sum(enrichedColsToAdd), nrow=nrow(enrichedGeneSet), data=0))
		colnames(emptyEnrichedToAdd)=paste0("enriched_in_", expectedThGeneSetsWithGeneSetAnalysis[enrichedColsToAdd])
	
		allColEnrichedGeneSet=cbind(enrichedGeneSet, emptyEnrichedToAdd)
		allColEnrichedGeneSet= allColEnrichedGeneSet[, c("GeneSet", paste0("enriched_in_", expectedThGeneSetsWithGeneSetAnalysis))]
	} else {
		allColEnrichedGeneSet= enrichedGeneSet
	}
	
	# add overlap_genes cols	
	tempAdd= data.frame(matrix(, ncol=ncol(allColEnrichedGeneSet)-1, nrow=nrow(allColEnrichedGeneSet)))
	colnames(tempAdd)=gsub("enriched", "overlap_genes", colnames(allColEnrichedGeneSet[2:ncol(allColEnrichedGeneSet)]))
	
	fullEnrichedGeneSet=cbind(allColEnrichedGeneSet, tempAdd)
	
	fullEnrichedGeneSet $totalThListsEnriched=rowSums(fullEnrichedGeneSet[, grepl("enriched_in", colnames(fullEnrichedGeneSet ))])
	
	fullEnrichedGeneSet = fullEnrichedGeneSet[order(fullEnrichedGeneSet $totalThListsEnriched, decreasing=TRUE),]


    ###
    ### identify gene sets that contain druggable genes
    ###
    fullEnrichedGeneSet$enrichedSetContainsDruggableGene=NA
    fullEnrichedGeneSet$anySetContainsDruggableGene=NA
    druggableGeneList= subset(geneByThGeneList, druggableGene )$gene
#    ThSetEnrichedCols=grep("enriched_in", colnames(fullEnrichedGeneSet), value=TRUE)
	ThGeneListEnrichedCols=grep("enriched_in", colnames(fullEnrichedGeneSet), value=TRUE)
	ThGeneListNamesFromCols=gsub("enriched_in_", "", ThGeneListEnrichedCols)
#     ThSetNamesFromCols=gsub("enriched_in_", "", ThSetEnrichedCols)
#    ThSet_inSetCols=grep("inSet_", colnames(geneByThGeneList), value=TRUE)
#    ThSetNamesFrom_inSetCols=gsub("inSet_", "", ThSet_inSetCols)
	ThGeneList_inSetCols =gsub("enriched_in_", "inSet_", ThGeneListEnrichedCols)
	ThSetEnrichedCols=grep("enriched_in", colnames(fullEnrichedGeneSet), value=TRUE)

	ThGeneListNamesFrom_inSetCols=gsub("inSet_", "", ThGeneList_inSetCols)


for (i in 1:nrow(fullEnrichedGeneSet)){
	# i=1
	thisGs=fullEnrichedGeneSet$GeneSet[i]
	enrichedInThGeneList= ThGeneListNamesFromCols [fullEnrichedGeneSet [i,ThGeneListEnrichedCols]==1]
	genesInThisPathway=subset(allPathwaysDF, GeneSet== thisGs)$Gene
	druggableGenesInThisPathway=subset(allPathwaysDF, GeneSet== thisGs & Gene %in% druggableGeneList)$Gene
	if (length(druggableGenesInThisPathway)!=0){
	  fullEnrichedGeneSet$anySetContainsDruggableGene[i]=TRUE
	  # test whether those genes are in a list with this gene set enriched
	  anyDruggableGenePresentInset= ThGeneListNamesFrom_inSetCols [colSums(subset(geneByThGeneList, gene %in% druggableGenesInThisPathway)[, ThGeneList_inSetCols])>0]
	  
	# for each GSEA gene set, get the list of genes overlapping with each ThLists 
	  theseGenesByThGeneList=subset(geneByThGeneList, gene %in% genesInThisPathway)
	  setsWithGenesInPathway= ThGeneList_inSetCols [colSums(theseGenesByThGeneList[, ThGeneList_inSetCols])>0] #ThGeneListNamesFrom_inSetCols
	  genesInPathways=unlist(lapply(setsWithGenesInPathway, function(x) paste(theseGenesByThGeneList$gene[as.logical(theseGenesByThGeneList[, x])], collapse=",")))
	  fullEnrichedGeneSet[i, gsub("inSet", "overlap_genes_in", setsWithGenesInPathway)]= genesInPathways

#    whichGenesPresentInSet=apply(

	  if (length(intersect(enrichedInThGeneList, anyDruggableGenePresentInset))>0){
	    fullEnrichedGeneSet$enrichedSetContainsDruggableGene[i]=TRUE
	    
	  }
	} else { # remove this pathway from enriched pathways since it's not druggable.
	  fullEnrichedGeneSet$anySetContainsDruggableGene[i]=FALSE
	  fullEnrichedGeneSet$enrichedSetContainsDruggableGene[i]=FALSE
	}
}
    fullEnrichedGeneSetWithDruggableThListGene=subset(fullEnrichedGeneSet, anySetContainsDruggableGene)[,grep("anySetContainsDruggableGene", colnames(fullEnrichedGeneSet), invert=TRUE, value=TRUE)]



    ### MAKE mega table

    system.time(genesInGeneSets<-merge(allPathwaysDF, geneByThGeneList, by.x="Gene", by.y="gene"))

    colnames(genesInGeneSets)=gsub("inSet_", "geneInSet_", colnames(genesInGeneSets))

    megaTable=merge(genesInGeneSets, fullEnrichedGeneSet, by="GeneSet")

    colnames(megaTable)=gsub("enriched_in_", "pathwayEnrichedInSet_", colnames(megaTable))

    megaTable$sumForRanking= rowSums(megaTable[,c("thListsContainingGene", "totalThListsEnriched")], na.rm=TRUE)

    megaTable= megaTable[order(megaTable$sumForRanking),]


    genesInThListsByGeneSets<-NULL
    for (i in 1:nrow(fullEnrichedGeneSet)){
        #i=1
        thisGs=enrichedGeneSet$GeneSet[i]
    #	dim(subset(megaTable, GeneSet==thisGs)); length(unique(subset(megaTable, GeneSet==thisGs)$Gene))
        genesInThListsByGeneSets=rbind(genesInThListsByGeneSets,data.frame(geneSet= thisGs, countInThLists=length(unique(subset(megaTable, GeneSet==thisGs)$Gene))))
    }

    ###
    ### Gene set overview and drug lists
    ###

    # add gene-specific information to gene set list
    fullEnrichedGeneSetWithDruggableThListGene$druggableGenesInThLists=unlist(lapply(fullEnrichedGeneSetWithDruggableThListGene $GeneSet, function(thisGs)  paste(unique(subset(genesInGeneSets,GeneSet %in% thisGs & druggableGene)$Gene), collapse=", ")))


    fullEnrichedGeneSetWithDruggableThListGene$allMemberGenesInThLists=unlist(lapply(fullEnrichedGeneSetWithDruggableThListGene $GeneSet, function(thisGs)  paste(unique(subset(genesInGeneSets,GeneSet %in% thisGs)$Gene), collapse=", ")))


    fullEnrichedGeneSetWithDruggableThListGene $drugs=unlist(lapply(fullEnrichedGeneSetWithDruggableThListGene $GeneSet, function(x) paste(unique(subset(dgDataInteractions, gene %in% subset(genesInGeneSets, GeneSet ==x)$Gene)$drug), collapse=", ")))

    fullEnrichedGeneSetWithDruggableThListGene$druggableGeneWithDrug =unlist(lapply(fullEnrichedGeneSetWithDruggableThListGene $GeneSet, function(x) paste(unique(subset(dgDataInteractions, gene %in% subset(genesInGeneSets, GeneSet ==x)$Gene)$rawInteraction), collapse=", ")))
    
    write.txt(fullEnrichedGeneSetWithDruggableThListGene, file=paste(gseaMetaNonEmpty$sampleName[1], "GeneSetAggregation.txt"))

    # druggable pathways
    dp=unique(megaTable[,c("GeneSet", "druggableGene")])
    druggableGeneSets=subset(dp, druggableGene)$GeneSet

    megaTable$druggablePathway= megaTable $GeneSet %in% druggableGeneSets

    dim(subset(megaTable, druggablePathway))

    ###
    ### Geneset TABLES
    ###

    commonGSEAcolsDF=unique(Reduce("rbind",lapply(gseaGeneSetListRaw, function(x) unique(x[,c("GeneSetName", "N_Genes_in_Gene_Set_(K)", "Description")]))))

    multiListGeneSets<-commonGSEAcolsDF


	## add support for values not present in gseaMetaNonEmpty
	#    for (i in 1:nrow(gseaMetaNonEmpty)){
    for (i in 1:nrow(gseaMeta)){
    		thisDataTag=gseaMeta$dataTag[i]
    		if (thisDataTag %in% names(gseaGeneSetListRaw)){
		    thisWide= gseaGeneSetListRaw[[thisDataTag]]
		    colnames(thisWide)[4:7]=paste0(thisDataTag, "_", colnames(thisWide)[4:7])
		    multiListGeneSets =merge(multiListGeneSets, thisWide[,c(1,4:7)], by="GeneSetName", all=TRUE)
		} else {
			genericColNames= colnames(gseaGeneSetListRaw[[1]][,4:7])
			emptyColsToAdd= data.frame(matrix(, ncol=length(genericColNames), nrow=nrow(multiListGeneSets)))
			colnames(emptyColsToAdd)=paste0(thisDataTag, "_", genericColNames)
			multiListGeneSets =cbind(multiListGeneSets, emptyColsToAdd)
	    }
	}


    # add contains druggable gene found in one gene list
    write.txt(multiListGeneSets, file=paste(gseaMetaNonEmpty $sampleName[1], "GeneSetDetailsPerList.txt"))


    ###
    ### create linh's overview
    ###
    colsToSkip=3

    pancanUp=subset(geneByThGeneList, inSet_pc_up==1)$gene
    pancanDown=subset(geneByThGeneList, inSet_pc_down==1)$gene
    pandiseaseUp=subset(geneByThGeneList, inSet_pd_up ==1)$gene
    pandiseaseDown=subset(geneByThGeneList, inSet_pd_down==1)$gene
    top5=subset(geneByThGeneList, inSet_top5 ==1)$gene
    pcdUp=intersect(pancanUp, pandiseaseUp)
    pcdDown=intersect(pancanDown, pandiseaseDown)

    pancanUp_DGIdb=unique(subset(dgDataInteractions, gene %in% pancanUp)$rawInteraction)
    pandiseaseUp_DGIdb=unique(subset(dgDataInteractions, gene %in% pandiseaseUp)$rawInteraction)
    pcdUp_DGIdb=unique(subset(dgDataInteractions, gene %in% pcdUp)$rawInteraction)
    top5_DGIdb=unique(subset(dgDataInteractions, gene %in% top5)$rawInteraction)

    colnames_pc=c("Pan-cancer Up", "Pan-cancer Down", "DGIDB", "TARGET", "GSEA", "k/K", "FDR")
    colnames_pcd= colnames_pc
    colnames_pcd[1:2]=c("PCD Up", "PCD Down")
    colnames_pd= colnames_pc
    colnames_pd[1:2]=c("Pan-disease Up", "Pan-disease Down")
    colnames_top5= colnames_pc[-2]
    colnames_top5[1]=c("top5")

    colnameList=list(colnames_pc, colnames_pcd, colnames_pd, colnames_top5)

    overviewlist=list(pancanUp, pancanDown, pancanUp_DGIdb, pcdUp, pcdDown, pcdUp_DGIdb, pandiseaseUp, pandiseaseDown, pandiseaseUp_DGIdb, top5, top5_DGIdb)

    overviewDF=data.frame(matrix(data="", nrow=1+max(sapply(overviewlist, length)), ncol=length(overviewlist)))

    for (i in 1:length(overviewlist)){
        thisText=overviewlist[[i]]
        if (length(thisText)>0)  overviewDF[2:(1+length(thisText)),i]= thisText
    }
    overviewDF[1,]= c(colnames_pc[1:3], colnames_pcd[1:3], colnames_pd[1:3], colnames_top5[1:2])

    overviewDF2=data.frame(matrix(data="", nrow=1+max(sapply(overviewlist, length)), ncol=length(unlist(colnameList))+3* colsToSkip))
    overviewDF2[1,]= c(colnames_pc, rep("", colsToSkip), colnames_pcd, rep("", colsToSkip), colnames_pd, rep("", colsToSkip), colnames_top5)

    overviewDF2[,1:3]= overviewDF[,1:3] # pc
    overviewDF2[,11:13]= overviewDF[,4:6] # pd
    overviewDF2[,21:23]= overviewDF[,7:9] # pcd
    overviewDF2[,31:32]= overviewDF[,10:11] # top5


    ###
    ### prioritized druggable genes
    ###

    # broadest list

    prioritizedDruggableGenes=druggableGeneByThGeneList



    # candidate genes are:
    # upoutlier or top five percent (not down outlier) genes marked druggable by dgidb 

    # consider whether the list is (yes/no)
    # in pancan up outliers
    # in pandisease up outliers

    #
    # names of gene-containing pathways that are enriched in a treehouse gene set (like pancan up outliers or pandisease up outliers)
    #
    pathwayByEnrichedList=lapply(ThSetEnrichedCols, function(ei) unlist(lapply(prioritizedDruggableGenes$druggableGene, function(x) paste(subset(allPathwaysDF, Gene ==x & GeneSet %in% fullEnrichedGeneSetWithDruggableThListGene$GeneSet[fullEnrichedGeneSetWithDruggableThListGene[,ei]==1])$GeneSet, collapse=","))))

    # namesOfEnrichedPathways=data.frame(t(cbindList(pathwayByEnrichedList)))
    #namesOfEnrichedPathways =data.frame(t(Reduce("cbind", pathwayByEnrichedList)), row.names=NULL)
    namesOfEnrichedPathways=data.frame(Reduce("cbind", pathwayByEnrichedList))

    colnames(namesOfEnrichedPathways)= ThSetEnrichedCols

    #
    # count of gene-containing pathways that are enriched in a treehouse gene set (like pancan up outliers or pandisease up outliers)
    #
    pathwayCountByEnrichedList=lapply(ThSetEnrichedCols, function(ei) unlist(lapply(prioritizedDruggableGenes$druggableGene, function(x) length(subset(allPathwaysDF, Gene ==x & GeneSet %in% fullEnrichedGeneSetWithDruggableThListGene$GeneSet[fullEnrichedGeneSetWithDruggableThListGene[,ei]==1])$GeneSet))))

    countsOfEnrichedPathways=data.frame(Reduce("cbind", pathwayCountByEnrichedList))
    colnames(countsOfEnrichedPathways)=paste0("count_of_", ThSetEnrichedCols)

    prioritizedDruggableGenesWithPathwayAndCount=cbind(prioritizedDruggableGenes, cbind(namesOfEnrichedPathways, countsOfEnrichedPathways)[,c(rbind(ThSetEnrichedCols, paste0("count_of_", ThSetEnrichedCols)))]) 


    # all enriched pathways 
    prioritizedDruggableGenesCorrespondingPathways =lapply(prioritizedDruggableGenes$druggableGene, function(x) subset(allPathwaysDF, Gene ==x & GeneSet %in% fullEnrichedGeneSetWithDruggableThListGene$GeneSet))

    prioritizedDruggableGenesWithPathwayAndCount$allEnrichedPathways=unlist(lapply(prioritizedDruggableGenesCorrespondingPathways, function(x) paste(x$GeneSet, collapse=",")))

    prioritizedDruggableGenesWithPathwayAndCount $countOfAllEnrichedPathways=unlist(lapply(prioritizedDruggableGenesCorrespondingPathways, function(x) length(x$GeneSet)))


    # is in a pathways that is enriched in the same treehouse gene sets the gene is in (like pancan up outliers or pandisease up outliers)



    # a drug corresponding to the gene has been recommended before

    # a drug corresponding to the gene has been recommended for a sample with one of the diseases identified by tumor map

    # is gene an up-pancan outlier because of tissue signal? e.g. is it a btk/flt3 in AML, in brain tumor, NGF -- # (gain a point -- not explained by tissue specificity
    # e.g. 129


    # sheet format:
    # add to druggable genes
    # drug-containing pathways that are enriched only in pancancer 
    # drug-containing pathways that are enriched only in pandisease
    # drug-containing pathways that are enriched in both pandisease and pancancer



    ###
    ### write amalgam excel workbook
    ###

    xlF=paste0(sampleID, "_gsea_dgidb_output.xlsx")

    xlsxInputList=list(overview=overviewDF2, GeneSetDetailsPerList= multiListGeneSets, GeneSetAggregation= fullEnrichedGeneSetWithDruggableThListGene, druggableGeneAggregation= prioritizedDruggableGenesWithPathwayAndCount, allGeneAggregation= geneByThGeneList)
    names(xlsxInputList)[1]=paste(sampleID, "overview")
    # write.xlsx(xlsxInputList,xlF, colNames=c(F, rep(T, length(xlsxInputList)-1)))

     write.txt(prioritizedDruggableGenesWithPathwayAndCount, file=paste(dgMeta$sampleName[1], "druggableGeneAggregation.txt")) ## moved from before addition of pathway counts

    ## # 
    # # is it in a treehouse gene list (pcup or pcdown) an enriched pathways

    # # is the drug target in that's enriched enriched in a pathways in gene list

    # # in the druggable gene list, find drugs that are in both pdUp and pcUp
    # druggableGeneByThGeneList

    # # if it's a drug we've recommended before that's good

    # # if it's a drug we've recommended before

    # # a drug corresponding to the gene has been recommended before

    # # a drug corresponding to the gene has been recommended for a sample with this disease before




    # library(openxlsx)
    # xlF=paste0(sampleID, "_gsea_dgidb_output.xlsx")

    # xlsxInputList=list(overview=overviewDF2, GeneSetDetailsPerList= multiListGeneSets, GeneSetAggregation= enrichedGeneSetWithDruggableThListGene, druggableGeneAggregation= druggableGeneByThGeneList, allGeneAggregation= geneByThGeneList)
    # names(xlsxInputList)[1]=paste(sampleID, "overview")
    # write.xlsx(xlsxInputList,xlF, colNames=c(F, rep(T, length(xlsxInputList)-1)))


}

#### Main ####

# dataDir="~/Documents/Dropbox/ucsc/projects/combine_pathways_drug_info/feb_27_ranked_drugs/TH01_0119_S01"
# pathwayFileLoc="~/Documents/Dropbox/ucsc/projects/gitCode/analysis-methods/union gsea and dgidb results"
# sampleID="TH01_0119_S01"


args<-commandArgs(TRUE)

sample.id<-args[1]
#sample.ids<-unlist(strsplit(args[1], " ", fixed=TRUE))
data.base.dir<-args[2]
pathway.file<-args[3]

combine.gsea.dgidb(sample.id, data.base.dir, pathway.file)

# combine.gsea.dgidb("TH03_0005_S01", "~/downloads/TH03_0005_S01", "~/downloads/msigdb.v5.2.symbols.gmt")
# combine.gsea.dgidb("C021_0006_001666", "~/downloads/TERTIARY_FAIL.C021_0006_001666", "~/downloads/msigdb.v5.2.symbols.gmt")

# sampleID="C021_0006_001666"; dataDir="~/downloads/TERTIARY_FAIL.C021_0006_001666";  pathwayFile="~/downloads/msigdb.v5.2.symbols.gmt"

# sampleID="TH03_0005_S01"; dataDir= "~/downloads/TH03_0005_S01";  pathwayFile="~/downloads/msigdb.v5.2.symbols.gmt"

# # for(sample in sample.ids){
# #     # Set data dir to the sample's data dir, and run the script
# #     data.dir<-paste(data.base.dir, sample, sep="/")
# #     combine.gsea.dgidb(sample, data.dir, pathway.file)
# # }