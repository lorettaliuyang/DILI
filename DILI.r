
generate_drug_type <- function()
{
    library("XML")
	drug_fda = read.table("Data/summary_04Q1_15Q3/DRUG_REAC_PT.txt", stringsAsFactors=F, sep="\t", quote="") ## get all drug in FAERS
	drug_fda = unique(drug_fda[,1])
	full_drugbank = xmlParse("E:/Public_database/DrugBank/full database.xml")  ## parse drugbank .xml file
	top = xmlRoot(full_drugbank)  ## find root node
	drug_list = c()
	type_list = c()
	for(i in 1:xmlSize(top))  ## xmlSize is to find total number of nodes
	{
	    dname = xmlValue(top[[i]][["name"]])  ## find the node "name" attribute
		dtype = as.vector(xmlAttrs(top[[i]])["type"])   ## find type
		drug_list = c(drug_list, dname)
		type_list = c(type_list, dtype)
	}
	drug_type_drugbank = cbind(tolower(drug_list), type_list)  ## make sure tolower drug name
	type_fda = drug_type_drugbank[match(drug_fda, drug_type_drugbank[,1]),2]
	drug_type_fda = cbind(drug_fda, type_fda)
	write.table(drug_type_fda, "Predictive_Model/model_v2/drug_type.txt", row.names=F, col.names=F, sep="\t", quote=F) 
}


### calculate reporting ratio and p-value ###
### a = Reports for drug of interest AND Reports for Event of Interest
### b = Reports for drug of interest AND Reports for All Other Events
### c = Reports for all other drugs AND Reports for Event of Interest
### d = Reports for all other drugs AND Reports for All Other Events
### N = a+b+c+d
### reporting_ratio = (a*N)/((a+c)*(a+b))
preprocess <- function()
{
    library("XML")
	library(preprocessCore)
	adr_interest = c("hepatobiliary disorders")  ## set up the phenotype of interest
	adr_interest_type = "SOC"  ## set up the hierarchy level of adr of interest
	drug_type_fda = read.table("Predictive_Model/model_v2/drug_type.txt", stringsAsFactors=F, sep="\t", quote="")
	
	drug_reac_pt =  read.table("Data/summary_04Q1_15Q3/DRUG_REAC_PT.txt", stringsAsFactors=F, sep="\t", quote = "")
	meddra = read.table("Data/ADR_hierarchy_processed.txt", sep="\t", quote="", stringsAsFactors=F, header=T, check.names=F)  ## read MedDRA hierarchy
	### process each level separately
	logp = c()
	RR = c()
	if(adr_interest_type=="PT")
	{
	    d_hie = unique(drug_reac_pt[which(adr_interest==drug_reac_pt[,2]),1])  ## get the drugs that can have adr matched to adr_interest
		f_hie = drug_reac_pt[ intersect(which(drug_reac_pt$V1 %in% d_hie), which(drug_reac_pt$V2==adr_interest)),4]
	    for(i in 1:length(d_hie))
		{
			x1 = f_hie[i]
			temp2 = which(d_hie[i]==drug_reac_pt$V1)
			x2 = sum(drug_reac_pt[temp2,4])-x1
			temp3 = which(adr_interest==drug_reac_pt$V2)
		    x3 = sum(drug_reac_pt[temp3,4])-x1
		    x4 = sum(drug_reac_pt[,4])-x1-x2-x3
			m = matrix(c(x1,x2,x3,x4),2,2)
		    pvalue = fisher.test(m, alternative="greater")$p.value
		    if(pvalue==0) logp = c(logp, -log10(1E-323)) else logp = c(logp, -log10(pvalue)) ###1E-323 is the smallest value R can handle  
			N=x1+x2+x3+x4
			RR = c(RR, (as.double(x1)*N)/(as.double(x1+x3)*(x1+x2)))  ## use as.double to transfer integer into double, to avoid interger overflow
						
		}
		d_not_hie = setdiff((unique(drug_reac_pt[,1])), d_hie)  ## get drugs that do not have hierarchy matched
		f_not_hie = rep(0, len=length(d_not_hie))
		for(i in 1:length(d_not_hie))
		{
			logp = c(logp, 0)
			RR = c(RR,0)
		}
		drug_matched = c(d_hie, d_not_hie)   ## get all drugs
		fre_matched = c(f_hie, f_not_hie)    ## get all frequencies
	}else if(adr_interest_type=="HLT")
	      {
		      drug_reac_hlt = drug_reac_pt
			  drug_reac_hlt[,2] = meddra[match(drug_reac_hlt[,2], meddra$"Term provided"),4]  ## map to hirerchy
			  drug_reac_hlt = aggregate(V4 ~ V1+V2+V3, data=drug_reac_hlt, sum)  ## sum the 4th column of rows if they have the same 1st, 2nd and 3rd column
			  d_hie = unique(drug_reac_hlt[which(adr_interest==drug_reac_hlt[,2]),1])  ## get the drugs that can have adr matched to adr_interest
		      f_hie = drug_reac_hlt[ intersect(which(drug_reac_hlt$V1 %in% d_hie), which(drug_reac_hlt$V2==adr_interest)),4]
			  for(i in 1:length(d_hie))
		      {
			      x1 = f_hie[i]
			      temp2 = which(d_hie[i]==drug_reac_pt$V1)  ## need drug_reac_pt instead of drug_reac_soc, since drug_reac_pt include all ADRs, and drug_reac_soc do not include un-matched ones
			      x2 = sum(drug_reac_pt[temp2,4])-x1
			      temp3 = which(adr_interest==drug_reac_hlt$V2)
		          x3 = sum(drug_reac_hlt[temp3,4])-x1
		          x4 = sum(drug_reac_pt[,4])-x1-x2-x3
			      m = matrix(c(x1,x2,x3,x4),2,2)
		          pvalue = fisher.test(m, alternative="greater")$p.value
		          if(pvalue==0) logp = c(logp, -log10(1E-323)) else logp = c(logp, -log10(pvalue)) ###1E-323 is the smallest value R can handle  
			      N=x1+x2+x3+x4
			      RR = c(RR, (as.double(x1)*N)/(as.double(x1+x3)*(x1+x2)))		
		      }
			  d_not_hie = setdiff((unique(drug_reac_hlt[,1])), d_hie)  ## get drugs that do not have hierarchy matched
			  f_not_hie = rep(0, len=length(d_not_hie))
			  for(i in 1:length(d_not_hie))
		      {
			      logp = c(logp, 0)
			      RR = c(RR,0)
		      }
			  drug_matched = c(d_hie, d_not_hie)   ## get all drugs
			  fre_matched = c(f_hie, f_not_hie)    ## get all frequencies			  
		  }else if(adr_interest_type=="HLGT")
		        {
		            drug_reac_hlgt = drug_reac_pt
			        drug_reac_hlgt[,2] = meddra[match(drug_reac_hlgt[,2], meddra$"Term provided"),5]
			        drug_reac_hlgt = aggregate(V4 ~ V1+V2+V3, data=drug_reac_hlgt, sum)
					d_hie = unique(drug_reac_hlgt[which(adr_interest==drug_reac_hlgt[,2]),1])  ## get the drugs that can have adr matched to adr_interest
		            f_hie = drug_reac_hlgt[ intersect(which(drug_reac_hlgt$V1 %in% d_hie), which(drug_reac_hlgt$V2==adr_interest)),4]
		            for(i in 1:length(d_hie))
		            {
			            x1 = f_hie[i]
			            temp2 = which(d_hie[i]==drug_reac_pt$V1)  ## need drug_reac_pt instead of drug_reac_soc, since drug_reac_pt include all ADRs, and drug_reac_soc do not include un-matched ones
			            x2 = sum(drug_reac_pt[temp2,4])-x1
			            temp3 = which(adr_interest==drug_reac_hlgt$V2)
		                x3 = sum(drug_reac_hlgt[temp3,4])-x1
		                x4 = sum(drug_reac_pt[,4])-x1-x2-x3
			            m = matrix(c(x1,x2,x3,x4),2,2)
		                pvalue = fisher.test(m, alternative="greater")$p.value
		                if(pvalue==0) logp = c(logp, -log10(1E-323)) else logp = c(logp, -log10(pvalue)) ###1E-323 is the smallest value R can handle  
			            N=x1+x2+x3+x4
			            RR = c(RR, (as.double(x1)*N)/(as.double(x1+x3)*(x1+x2)))		
		            }
					d_not_hie = setdiff((unique(drug_reac_hlgt[,1])), d_hie)  ## get drugs that do not have hierarchy matched
					f_not_hie = rep(0, len=length(d_not_hie))
					for(i in 1:length(d_not_hie))
		            {
			            logp = c(logp, 0)
			            RR = c(RR,0)
		            }
					drug_matched = c(d_hie, d_not_hie)   ## get all drugs
					fre_matched = c(f_hie, f_not_hie)    ## get all frequencies
				}else
				{
				    drug_reac_soc = drug_reac_pt
			        drug_reac_soc[,2] = meddra[match(drug_reac_soc[,2], meddra$"Term provided"),6]
			        drug_reac_soc = aggregate(V4 ~ V1+V2+V3, data=drug_reac_soc, sum)
					d_hie = unique(drug_reac_soc[which(adr_interest==drug_reac_soc[,2]),1])  ## get the drugs that can have adr matched to adr_interest
		            f_hie = drug_reac_soc[ intersect(which(drug_reac_soc$V1 %in% d_hie), which(drug_reac_soc$V2==adr_interest)),4]
					for(i in 1:length(d_hie))
					{
					    x1 = f_hie[i]
						temp2 = which(d_hie[i]==drug_reac_pt$V1)   ## need drug_reac_pt instead of drug_reac_soc, since drug_reac_pt include all ADRs, and drug_reac_soc do not include un-matched ones
						x2 = sum(drug_reac_pt[temp2,4])-x1
						temp3 = which(adr_interest==drug_reac_soc$V2)
		                x3 = sum(drug_reac_soc[temp3,4])-x1
		                x4 = sum(drug_reac_pt[,4])-x1-x2-x3
						m = matrix(c(x1,x2,x3,x4),2,2)
		                pvalue = fisher.test(m, alternative="greater")$p.value
						if(pvalue==0) logp = c(logp, -log10(1E-323)) else logp = c(logp, -log10(pvalue)) ###1E-323 is the smallest value R can handle    
						N=x1+x2+x3+x4
						RR = c(RR, (as.double(x1)*N)/(as.double(x1+x3)*(x1+x2)))
						
					}
					
					d_not_hie = setdiff((unique(drug_reac_soc[,1])), d_hie)  ## get drugs that do not have hierarchy matched
					f_not_hie = rep(0, len=length(d_not_hie))
					for(i in 1:length(d_not_hie))
					{
					    logp = c(logp, 0)
						RR = c(RR,0)
					}
					drug_matched = c(d_hie, d_not_hie)   ## get all drugs
					fre_matched = c(f_hie, f_not_hie)    ## get all frequencies
				}
	
	drug_modeling <- intersect(drug_type_fda[which(drug_type_fda$V2=="small molecule"),1], drug_matched)  ## modeling drug are small molecule AND have adr reported in FDA
	drug_modeling <- sort(drug_modeling)  ## make sure the drug name are ordered
	fre_modeling <- fre_matched[match(drug_modeling, drug_matched)]	 ## get corresponding fre for modeling drugs
	logp_modeling <- logp[match(drug_modeling, drug_matched)]
	RR_modeling <- RR[match(drug_modeling, drug_matched)]
	
	## get the patent year information from drugbank
	full_drugbank = xmlParse("E:/Public_database/DrugBank/full database.xml")  ## parse drugbank .xml file
	top = xmlRoot(full_drugbank)  ## find root node
	drug_list = c()
	year_list = c()
	for(i in 1:xmlSize(top))  ## xmlSize is to find total number of nodes
	{
	    dname = tolower(xmlValue(top[[i]][["name"]]))  ## find the node "name" attribute
		allyear = c()
		if((xmlValue(top[[i]][["patents"]]))!="")
		{
		    for(j in 1:xmlSize(top[[i]][["patents"]]))   ## find how many patent nodes in patents
			{
			    dates = xmlValue(top[[i]][["patents"]][[j]][["approved"]])   ## find approved date
				allyear = c(allyear, strsplit(dates, "-")[[1]][1])   ## find approved year for different patnet node
			}
			year = min(allyear)   ## find the earliest approved year
		}else
		{
		    year = "N"
		}
		drug_list = c(drug_list, dname)
		year_list = c(year_list, year)
	}
	year_modeling <- year_list[match(drug_modeling, drug_list)]
	
    ## merge all information together	
	drug_fp =  read.csv("drug_fp/descriptors.csv", header=TRUE, stringsAsFactors=F, sep=",", quote="", check.names=F)   ## fp file generated by MOE
	drugbank <- read.table("E:/Public_database/DrugBank/id2name.txt", stringsAsFactors=F, sep="\t", header=T, quote = "") 
	d_id = drugbank[match(drug_modeling, tolower(drugbank$Name)), "ID"]
	drug_modeling_feature = drug_fp[match(d_id, drug_fp$ID),3:820]   ## get descriptor info column 3 - column 820
	drug_modeling_feature = cbind(drug_modeling, fre_modeling, logp_modeling, RR_modeling, year_modeling, drug_modeling_feature)   ## put the 1st and 2nd column as the drug name and frequency
	drug_modeling_feature = na.omit(drug_modeling_feature)   ## only get the non-NA rows
	write.table(drug_modeling_feature, "Predictive_Model/model_v2/drug_modeling_feature_2.txt", row.names=F, col.names=T, sep="\t", quote=F)
	
}
correlation_stat <- function()
{
    mfeature = read.table("Predictive_Model/model_v2/results_v2.4/topleft_3/drug_modeling_feature_normalized_2.txt",  sep="\t", quote="", stringsAsFactors=F, header=T, comment.char="", check.names=F)
	drug_modeling = mfeature$drug_modeling
	fre_modeling = mfeature$fre_modeling
	
	mlabel = rep(1, len=length(drug_modeling))  ## get the label 
	mlabel[which(fre_modeling==0)] = -1
	
	results = apply(as.matrix(mfeature[,6:ncol(mfeature)]), 2, function(x) cor.test(x,mlabel))
	feature = names(results)
	pvalue = as.vector(sapply(results, function(x) x$p.value))
	correlation = as.vector(sapply(results, function(x) x$estimate))
	write.table(cbind(feature, pvalue, correlation), "Predictive_Model/model_v2/results_v2.4/topleft_3/drug_correlation.txt", row.names=F, col.names=T, sep="\t", quote=F)
	
	feature_sig = feature[which(pvalue<=0.05)]
	pvalue_sig = pvalue[which(pvalue<=0.05)]
	correlation_sig = correlation[which(pvalue<=0.05)]
	write.table(cbind(feature_sig, pvalue_sig, correlation_sig), "Predictive_Model/model_v2/results_v2.4/topleft_3/drug_correlation_sig.txt", row.names=F, col.names=T, sep="\t", quote=F)
}

### this function works on all features first, rank them and get top 100 features, then work on these 100 features, rank them, filter out bottom 10, then work on 90 features, rank them filter out bottom 10, so on.
### in each test, also tune cost to get the best performance
### work on a fixed data set: all drugs with frequency equal to 0 will be negative, drugs with frequency>=20 will be positive, sapling without replacement to keep balanced
### changes compared to testing_3: 
### 1) normalize data to -1 and 1 before RFE and SVM, delete features that with Na values (those are the features that are with the same value to all drugs);
### 2) add drug name and fre in feature file
### 3) correct the feature name in model file
### 4) fixed the cost parameter in RFE
### 5) normalize the data based on pos and neg, not based on all drugs
training_testing <- function(x)
{
    library(sampling)
	library(e1071)
	library(OmicsMarkeR)
	library(pROC)
	library(clusterSim)
	current_year = 2016
	
	feature_data = read.table("Predictive_Model/model_v2/drug_modeling_feature_2.txt",  sep="\t", quote="", stringsAsFactors=F, header=T, comment.char="", check.names=F)
	drug_all = feature_data$drug_modeling
	fre_all = feature_data$fre_modeling
	logp_all = feature_data$logp_modeling
	RR_all = feature_data$RR_modeling
	year_all = feature_data$year_modeling
	
	## pos>=20, neg=0, sampling without replacement to keep balanced
	pos_index = intersect(intersect(which(fre_all>=15), which(logp_all>1.3)), which(RR_all>2))  ## pos is drugs that satisfy multiple criterial, logp of 1.3 equal to p of 0.05
	neg_index = intersect(which(fre_all==0), which(year_all<(current_year-10)))   ## neg is drugs keep 0 reports for at least 10 years
	all_index = c(pos_index, neg_index)
	drug_modeling = feature_data[all_index, "drug_modeling"]
	fre_modeling = feature_data[all_index, "fre_modeling"]
	logp_modeling = feature_data[all_index, "logp_modeling"]
	RR_modeling = feature_data[all_index, "RR_modeling"]
	year_modeling = feature_data[all_index, "year_modeling"]
	feature = feature_data[all_index,6:ncol(feature_data)]
	feature = data.Normalization(feature,type="n5")  ## need package of clusterSim, n5 normalize data between -1 and 1
	feature = feature[,colSums(is.na(feature))==0]   ##remove the columns that have at least one NA, colSums count the number of columns that meet certain criterial
	write.table(cbind(drug_modeling, fre_modeling, logp_modeling, RR_modeling, year_modeling, feature), "Predictive_Model/model_v2/results_v2.4/topleft_15_1.3_2/drug_modeling_feature_normalized.txt", row.names=F, col.names=T, sep="\t", quote=F)
	
	### keep pos and neg balanced
	testing_index = union(sample(pos_index, round(length(pos_index)*0.1), replace=F), sample(neg_index, round(length(pos_index)*0.1), replace=F))   ### round: half adjust. get testing drug 10% from all pos, and 10% from all neg
	drug_testing = feature_data[testing_index, "drug_modeling"]
	drug_testing = sort(drug_testing)
	fre_testing = fre_modeling[match(drug_testing, drug_modeling)]
	label_testing = rep("pos", len=length(drug_testing))
	label_testing[which(fre_testing==0)] = "neg"
	feature_testing = feature[match(drug_testing, drug_modeling),]	## get corresponding fre for testing drugs
	#feature_testing = cbind(feature_testing, label_testing)

	drug_training_temp = setdiff(drug_modeling, drug_testing)  ## training drug are small molecule AND have adr reported in adr_interest (drug_matched), but not testing drug. first find drugs that are not in testing
    ## pos in training will be intersect drugs in both of drug_modeling[pos] and temp; neg in training will random selected from intersect of drug_modeling[neg] and temp, with the size equal to pos
	drug_training = c( intersect(drug_training_temp,drug_all[pos_index]), sample(intersect(drug_training_temp,drug_all[neg_index]), length(intersect(drug_training_temp,drug_all[pos_index])), replace=F))  ## keep pos and neg balanced 
	drug_training = sort(drug_training)  ## make sure the drug name are ordered
	fre_training = fre_modeling[match(drug_training, drug_modeling)]
	label_training = rep("pos", len=length(drug_training))
	label_training[which(fre_training==0)] = "neg"
	feature_training = feature[match(drug_training, drug_modeling),]	## get corresponding fre for training drugs
	#feature_training = cbind(feature_training, label_training)
	
	filename1 = paste("Predictive_Model/model_v2/results_v2.4/topleft_15_1.3_2/test_fixed_RFE_", x, sep="")
	filename1 = paste(filename1, "drug_training_feature.txt", sep="/")
	filename2 = paste("Predictive_Model/model_v2/results_v2.4/topleft_15_1.3_2/test_fixed_RFE_", x, sep="")
	filename2 = paste(filename2, "drug_testing_feature.txt", sep="/")
	write.table(cbind(drug_training, fre_training, label_training, feature_training), filename1, row.names=F, col.names=F, sep="\t", quote=F)
	write.table(cbind(drug_testing, fre_testing, label_testing, feature_testing), filename2, row.names=F, col.names=F, sep="\t", quote=F)
	
	top_threshold = c( ncol(feature), 100, 90,80,70,60,50,40,30,20,15,10)
	cost_threshold = c(1, 2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100)
		
	auc_value = c()   ## iniatiate the auc column in auc file
	cost_seq = c()    ## iniatiate the cost column in auc file
	top_seq = c()   ## iniatiate the top column in auc file
	for(k in 1:length(cost_threshold))
	{
	    filename = paste("Predictive_Model/model_v2/results_v2.4/topleft_15_1.3_2/test_fixed_RFE_", x, sep="")
		filename = paste(filename, "ROC_curve_", sep="/")
		filename = paste(filename, cost_threshold[k], sep="")
		filename = paste(filename, ".pdf", sep="")
		pdf(file=filename, width=10, height=10)   ## for each cost threshold, draw one figure that include multiple plots, width and height define the figure size
		par(mfrow = c(3,4), mar=c(1,1,1,1))   ## multiple plots, mar defines the margin
		
		selected_name = colnames(feature)  ## initially use all features
		for(j in 1:length(top_threshold))
		{
		    svm.model = svm(feature_training[,selected_name], factor(label_training), cost=cost_threshold[k], kernel="linear", scale=F, cross=10, probability=TRUE)
		    svm.pred = predict(svm.model, feature_testing[,selected_name], decision.values = T, probability=T) 
			
			## write model summary
	        filename = paste("Predictive_Model/model_v2/results_v2.4/topleft_15_1.3_2/test_fixed_RFE_", x, sep="")
			filename = paste(filename, "results_training_model_svm_", sep="/")
			filename = paste(filename, cost_threshold[k], sep="")
			filename = paste(filename, top_threshold[j], sep="_")
			filename = paste(filename, ".txt", sep="")
	        write.table(capture.output(summary(svm.model)), filename, row.names=F,  col.names=F, quote=F)  ## write the best model 
	        
			## write selected feature
	        filename = paste("Predictive_Model/model_v2/results_v2.4/topleft_15_1.3_2/test_fixed_RFE_", x, sep="")
			filename = paste(filename, "results_selected_features_", sep="/")
			filename = paste(filename, cost_threshold[k], sep="")
			filename = paste(filename, top_threshold[j], sep="_")
			filename = paste(filename, ".txt", sep="")
			weight <- t(svm.model$coefs) %*% svm.model$SV   ## weight is a matrix: row=1, col=number of features
			write.table(cbind(selected_name, as.vector(weight)), filename, row.names=T,  col.names=F, quote=F, sep="\t")  ## write the best model, if directly output weight value, some of the feature names with special characters would be incorrect, so output feature name and weight value separately using cbind
			
            ## write the model before update the selected_name, because if update selected_name first, the selected_name and weight value would not be in the same length		
	        if(j!=length(top_threshold)) ## when j==length(top_threshold), no need to rank the feature any more
			{
			    ranked_index = svmrfeFeatureRanking(feature_training[,selected_name], factor(label_training), c=cost_threshold[k], perc.rem=10)   ## use the feature name instead of feature index to rfe, because rfe output index in new feature list instead of their index in the original data
		        ranked_name = selected_name[ranked_index]
				selected_name = ranked_name[1:top_threshold[j+1]]  ## update selected_name to the top ones selected
			}
			
	        ## get true label and predict label, calculate performance
	        label_assign = label_testing  ## get the true label
	        label_probability = as.vector(attr(svm.pred, "probabilities")[,"pos"]) ## get the predicting probability of being pos
			
			### write the roc curve
			roc_curve = roc(label_assign, label_probability, levels=c("neg", "pos"))  ## levels help to define control(neg in this case) and case (pos in this case) samples in the data   
			## auc.polygon: whether fill in the roc area; grid: setup grid in the plot with interval; col: the color of roc curve
			plot(roc_curve, print.thres=T, print.auc=T, col="blue", auc.polygon=T, grid=c(0.1,0.2), grid.col=c("green", "red"), auc.polygon.col="lightcyan")
			t1 = paste("cost=", cost_threshold[k], sep="")
			t2 = paste("top=", top_threshold[j], sep="")
			text_value = c(t1, t2)
			legend("topleft", legend=text_value, text.font=4)  ## text.font: an integer specifying the font style of the legend text; possible values are: normal(1), bold(2), italic(3), bold and italic(4)
            
            ### get the auc value			
			auc_value = c(auc_value, auc(roc_curve))   ## prepare the auc column in auc file
			cost_seq = c(cost_seq, cost_threshold[k])  ## prepare the cost column in auc file
			top_seq = c(top_seq, top_threshold[j])    ## prepare the top column in auc file
			
			## write the performance under best threshold, and also performance under all thresholds
			best_threshold = coords(roc_curve, "b", ret=c("threshold"), best.method=c("closest.topleft")) ## coords: output coordinate; "b": best threshold under the method of topleft; ret: output parameters
			if(length(best_threshold)>1)  ## if two points have exact the same sentivity and specificity (may be in reverse in the two situations), choose the point that has larger sensitivity
			{
			    temp = coords(roc_curve, best_threshold, ret=c("threshold", "sensitivity"), best.method=c("closest.topleft"))  ## output threshold and sensitivity under the situations of best_threshold
				best_threshold = temp["threshold", which(temp["sensitivity",]== max(temp["sensitivity",]))]  ## reassign best_threhold to the one with higher sensitivity
			}
			best_sensitivity = coords(roc_curve, best_threshold, ret=c("sensitivity"), best.method=c("closest.topleft"))
			best_specificity = coords(roc_curve, best_threshold, ret=c("specificity"), best.method=c("closest.topleft"))
			best_FPR = coords(roc_curve, best_threshold, ret=c("1-specificity"), best.method=c("closest.topleft"))
			best_FNR = coords(roc_curve, best_threshold, ret=c("1-sensitivity"), best.method=c("closest.topleft"))
			best_accuracy = coords(roc_curve, best_threshold, ret=c("accuracy"), best.method=c("closest.topleft"))
			best_PPV = coords(roc_curve, best_threshold, ret=c("ppv"), best.method=c("closest.topleft"))
			bestvalue = c(cost_threshold[k], top_threshold[j], best_threshold, best_sensitivity, best_specificity, best_FPR, best_FNR, best_accuracy, best_PPV)
			names(bestvalue) = c("cost_threshold", "top_threshold", "best_probability", "best_sensitivity", "best_specificity", "best_FPR", "best_FNR", "best_accuracy", "best_PPV")
			allvalues = coords(roc=roc_curve, roc_curve$threshold, input="threshold", ret=c("threshold","sensitivity", "specificity", "accuracy", "tn", "tp", "fn", "fp","npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv"))
            filename = paste("Predictive_Model/model_v2/results_v2.4/topleft_15_1.3_2/test_fixed_RFE_", x, sep="")
			filename = paste(filename, "results_testing_performance_svm_", sep="/")
			filename = paste(filename, cost_threshold[k], sep="")
	        filename = paste(filename, top_threshold[j], sep="_")
			filename = paste(filename, ".txt", sep="")
			write.table(bestvalue, filename, row.names=T, col.names=F, sep="\t", quote=F, append=T)  ## output training performance
	        write.table(allvalues, filename, row.names=T, col.names=F, sep="\t", quote=F, append=T)  ## output training performance
			
			### write the testing label in the situation of best threshold
	        label_assign = label_testing  ## get the true label
	        label_predict = rep("neg", length(label_assign))  ## create "neg" vector
	        label_predict[which(label_probability>=best_threshold)] = "pos"
	        mlabel = cbind(drug_testing, fre_testing, label_assign, label_predict, label_probability)
	        colnames(mlabel) = c("Drug", "Frequency", "Real", "Predict", "Probability")
	        filename = paste("Predictive_Model/model_v2/results_v2.4/topleft_15_1.3_2/test_fixed_RFE_", x, sep="")
			filename = paste(filename, "results_testing_label_svm_", sep="/")
			filename = paste(filename, cost_threshold[k], sep="")
	        filename = paste(filename, top_threshold[j], sep="_")
			filename = paste(filename, ".txt", sep="")
	        write.table(mlabel, filename, row.names=F, col.names=T, sep="\t", quote=F)
		}
	    dev.off()
	}
	mroc = cbind(cost_seq, top_seq, auc_value)
	colnames(mroc) = c("cost", "top", "auc")
	filename = paste("Predictive_Model/model_v2/results_v2.4/topleft_15_1.3_2/test_fixed_RFE_", x, sep="")
	filename = paste(filename, "results_auc.txt", sep="/")
	write.table(mroc, filename, row.names=F, col.names=T, sep="\t", quote=F)
}

### read through all auc.txt file, get the average value of different data sets under the same cost value and top value, find the best cost and top value
summarize_auc <- function()
{
	auc_value = c()
	for(i in 1:times)
	{
	    filename = paste("Predictive_Model/model_v2/results_v2.4/topleft_15_1.3_2/test_fixed_RFE_", i, sep="")
		filename = paste(filename, "results_auc.txt", sep="/")
		auc_file = read.table(filename,  sep="\t", quote="", stringsAsFactors=F, header=T, check.names=F)
		auc_value = cbind(auc_value, auc_file[,3])  ## include each auc value
		
	}
	ave_auc = rowMeans(auc_value)   ## calculate the element-wise row mean of a matrix
	cost = auc_file[,1]
	top = auc_file[,2]
	mdata = cbind(cost,top,ave_auc)
	colnames(mdata) = c("cost", "top", "auc")
	write.table(mdata, "Predictive_Model/model_v2/results_v2.4/topleft_15_1.3_2/average_auc.txt", row.names=F, col.names=T, sep="\t", quote=F)
	
	## find the best cost, top, auc
	best_data = mdata[order(-mdata[,"auc"]),]   ## order mdata based on column of auc
	best_cost <<- best_data[1, "cost"]
	best_top <<- best_data[1, "top"]
	best_auc <<- best_data[1, "auc"]
}

### read through all performance file under the condition of best cost value and best top value, get the average sensitivity, specificity, probability threshold for different data sets
### write all results into model_overall_performance.txt
summarize_performance <- function()
{
	sensitivity = c()
	specificity = c()
	FPR = c()
	FNR = c()
	accuracy = c()
	PPV = c()
	probability = c()
	for(i in 1:times)
	{
	    filename = paste("Predictive_Model/model_v2/results_v2.4/topleft_15_1.3_2/test_fixed_RFE_", i, sep="")
		filename = paste(filename, "results_testing_performance_svm_", sep="/")
		filename = paste(filename, best_cost, sep="")
		filename = paste(filename, best_top, sep="_")
		filename = paste(filename, ".txt", sep="")
		perf_file = read.table(filename,  sep="\t", quote="", stringsAsFactors=F, nrows=9)   ## only read the first 9 rows
		sensitivity = c(sensitivity, perf_file[which(perf_file$V1=="best_sensitivity"), 2])
		specificity = c(specificity, perf_file[which(perf_file$V1=="best_specificity"), 2])
		FPR = c(FPR, perf_file[which(perf_file$V1=="best_FPR"), 2])
		FNR = c(FNR, perf_file[which(perf_file$V1=="best_FNR"), 2])
		accuracy = c(accuracy, perf_file[which(perf_file$V1=="best_accuracy"), 2])
		PPV = c(PPV, perf_file[which(perf_file$V1=="best_PPV"), 2])
		probability = c(probability, perf_file[which(perf_file$V1=="best_probability"), 2])
	
	}	
    sensitivity = mean(sensitivity)	
	specificity = mean(specificity)
	FPR = mean(FPR)
	FNR = mean(FNR)
	accuracy = mean(accuracy)
	PPV = mean(PPV)
	probability = mean(probability)

	model_perf = c(best_cost, best_top, sensitivity, specificity, FPR, FNR, accuracy, PPV, probability)
	names(model_perf) = c("cost", "top", "sensitivity", "specificity", "FPR", "FNR", "accuracy", "PPV", "probability")
    write.table(model_perf, "Predictive_Model/model_v2/results_v2.4/topleft_15_1.3_2/model_overall_performance.txt", row.names=T, col.names=F, sep="\t", quote=F)
} 

summarize_feature <- function()
{
    library(plyr)
	feature_all = c()
	for(i in 1:times)
	{
	    filename = paste("Predictive_Model/model_v2/results_v2.4/topleft_15_1.3_2/test_fixed_RFE_", i, sep="")
		filename = paste(filename, "results_selected_features_", sep="/")
		filename = paste(filename, best_cost, sep="")
		filename = paste(filename, best_top, sep="_")
		filename = paste(filename, ".txt", sep="")
		feature_file = read.table(filename, sep="\t", quote="", stringsAsFactors=F, comment.char="", check.names=F)  
		feature_all = c(feature_all, feature_file[,2])
	}
	feature_fre = count(feature_all) 
	feature_fre = feature_fre[order(-feature_fre[,"freq"]),]
	colnames(feature_fre) = c("feature", "frequency")
	write.table(feature_fre[,1:2], "Predictive_Model/model_v2/results_v2.4/topleft_15_1.3_2/selected_features.txt", row.names=F, col.names=T, sep="\t", quote=F)
}


### genrate subset using selected features, all pos and all neg
subset_data <- function()
{
    current_year = 2016
	feature_data = read.table("Predictive_Model/model_v2/results_v2.4/topleft_3/drug_modeling_feature_normalized_2.txt",  sep="\t", quote="", stringsAsFactors=F, header=T, comment.char="", check.names=F)
	drug_modeling = feature_data$drug_modeling
	fre_modeling = feature_data$fre_modeling
	logp_modeling = feature_data$logp_modeling
	RR_modeling = feature_data$RR_modeling
	year_modeling = feature_data$year_modeling
	feature = feature_data[,6:ncol(feature_data)]   ## when read normalized data, the number of features decreased since some features have been eleminated during pre-process
	
	## define the rows to be selected
	pos_index = intersect(intersect(which(fre_modeling>=20), which(logp_modeling>2)), which(RR_modeling>2))  ## pos is drugs that satisfy multiple criterial, logp of 1.3 equal to p of 0.05
	neg_index = intersect(which(fre_modeling==0), which(year_modeling<(current_year-10)))   ## neg is drugs keep 0 reports for at least 10 years
	## define the columns to be selected
	selected_feature = c("drug_modeling", "fre_modeling", "logp_modeling", "RR_modeling", "year_modeling", "a_ICM", "key62", "key335", "key297", "SMR_VSA7", "key25", "key346", "key15", "SSF043", "key137")
	sub_feature_data = feature_data[c(pos_index,neg_index),selected_feature]
	write.table(sub_feature_data, "Predictive_Model/model_v2/results_v2.4/topleft_3/subset_drug_modeling_feature_normalized_2.txt", row.names=F, col.names=T, sep="\t", quote=F)
	
}

heatmap_draw <- function()
{
    mdata = read.table("Predictive_Model/model_v2/drug_modeling_feature_normalized.txt",  sep="\t", quote="", stringsAsFactors=F, header=T, comment.char="", check.names=F)
	selected_feature = read.table("E:/Drug_ADR_Networks/Predictive_Model/model_v2/results_v2.4/selected_features.txt", stringsAsFactors=F, header=T, sep="\t")
	selected_feature = selected_feature[1:30,1]
	mdata_part = mdata[, c("drug_modeling", "fre_modeling", selected_feature)]
    write.table(mdata_part, "Predictive_Model/model_v2/results_v2.4/drug_modeling_selected_features.txt", row.names=F, col.names=T, sep="\t", quote=F)
	
	heatdata = read.table("Predictive_Model/model_v2/results_v2.4/drug_modeling_selected_features_0_5.txt", sep="\t", quote="", stringsAsFactors=F, header=T, comment.char="", check.names=F)
	weight = heatdata[,3:32]  ## weight is all feature columns
	weight = as.matrix(weight)
	## Rowv and Colv indicate wether cluster the row and column, default is clustering, labRow is the row name,
	mheatmap <- heatmap(weight, Rowv=NA, Colv=NA, col = heat.colors(100), margins=c(5,10), scale="column", labRow=heatdata[,1])
	mheatmap <- heatmap(weight, Rowv=NA, col = heat.colors(100), margins=c(5,10), scale="column", labRow=heatdata[,1])
}

### compare different normalization methods: try different normalization methods, compare correlation between original data and normalized data based on some selected features
compare_norm <- function()
{
    library(sampling)
	library(e1071)
	library(OmicsMarkeR)
	library(clusterSim)
	current_year = 2016
	
	feature_data = read.table("Predictive_Model/model_v2/drug_modeling_feature_2.txt",  sep="\t", quote="", stringsAsFactors=F, header=T, comment.char="", check.names=F)
	drug_all = feature_data$drug_modeling
	fre_all = feature_data$fre_modeling
	logp_all = feature_data$logp_modeling
	RR_all = feature_data$RR_modeling
	year_all = feature_data$year_modeling
	
	## pos>=20, neg=0, sampling without replacement to keep balanced
	pos_index = intersect(intersect(which(fre_all>=20), which(logp_all>2)), which(RR_all>2))  ## pos is drugs that satisfy multiple criterial, logp of 1.3 equal to p of 0.05
	neg_index = intersect(which(fre_all==0), which(year_all<(current_year-10)))   ## neg is drugs keep 0 reports for at least 10 years
	all_index = c(pos_index, neg_index)
	drug_modeling = feature_data[all_index, "drug_modeling"]
	fre_modeling = feature_data[all_index, "fre_modeling"]
	logp_modeling = feature_data[all_index, "logp_modeling"]
	RR_modeling = feature_data[all_index, "RR_modeling"]
	year_modeling = feature_data[all_index, "year_modeling"]
	selected_feature = c("a_ICM", "SSF083", "key335", "SS036", "key62", "key307", "SSF043", "key25", "key137", "SSF068", "key15", "SMR_VSA7", "SS023", "key297", "key33", "key339", "key318", "kS_aaO", "key252", "SSF002", "SS019", "key13", "key200", "key108", "key28", "SS004", "key155", "ADF101", "SlogP_VSA9", "SSF015")
	feature = feature_data[all_index,selected_feature]
	feature1 = data.Normalization(feature,type="n1") 
	cor1 = diag(cor(feature, feature1))
	feature2 = data.Normalization(feature,type="n2") 
	cor2 = diag(cor(feature, feature2))
	feature3 = data.Normalization(feature,type="n3") 
	cor3 = diag(cor(feature, feature3))
	feature4 = data.Normalization(feature,type="n4") 
	cor4 = diag(cor(feature, feature4))
	feature5 = data.Normalization(feature,type="n5") 
	cor5 = diag(cor(feature, feature5))
	feature6 = data.Normalization(feature,type="n6") 
	cor6 = diag(cor(feature, feature6))
	feature12 = data.Normalization(feature,type="n12") 
	cor12 = diag(cor(feature, feature12))
	feature13 = data.Normalization(feature,type="n13") 
	cor13 = diag(cor(feature, feature13))
	mcomp = rbind(cor1,cor2,cor3,cor4,cor5,cor6,cor12,cor13)
	colnames(mcomp) = selected_feature
	write.table(mcomp, "Predictive_Model/model_v2/normalization_methods_compare.txt", sep="\t", col.names=T, row.names=T, quote=F)
}


###### use the selected best cost, best top, best features to train all pos and neg, and get one final model, use this final model to predict any independent data ######
independent_testing <- function()
{
    library(e1071)
	library(pROC)
	library(clusterSim)
	current_year = 2016
	### get all drugs and features
	feature_data = read.table("Predictive_Model/model_v2/drug_modeling_feature_2.txt",  sep="\t", quote="", stringsAsFactors=F, header=T, comment.char="", check.names=F)
	drug_modeling = feature_data$drug_modeling
	fre_modeling = feature_data$fre_modeling
	logp_modeling = feature_data$logp_modeling
	RR_modeling = feature_data$RR_modeling
	year_modeling = feature_data$year_modeling
	
	pos_index = intersect(intersect(which(fre_modeling>=15), which(logp_modeling>1.3)), which(RR_modeling>2))  ## pos is drugs that satisfy multiple criterial, logp of 1.3 equal to p of 0.05
	neg_index = intersect(which(fre_modeling==0), which(year_modeling<(current_year-10)))   ## neg is drugs keep 0 reports for at least 10 years
	
	## load the best performance
	performance = read.table("Predictive_Model/model_v2/results_v2.4/topleft_15_1.3_2/model_overall_performance.txt", sep="\t", quote="", stringsAsFactors=F)
	best_cost = performance[which(performance$V1=="cost"),2]
    best_top = performance[which(performance$V1=="top"),2]
	best_probability = performance[which(performance$V1=="probability"),2]
	selected_feature_all = read.table("Predictive_Model/model_v2/results_v2.4/topleft_15_1.3_2/selected_features.txt", sep="\t", header=T, quote="", stringsAsFactors=F, check.names=F, comment.char="")
	selected_feature = selected_feature_all[1:best_top,1]  ## could auto load top features, could also manually input the best top features
	##selected_feature = c("a_ICM", "SSF083", "key335", "SS036", "key62", "key307", "SSF043", "key25", "key137", "SSF068", "key15", "SMR_VSA7", "SS023", "key297", "key33", "key339", "key318", "kS_aaO", "key252", "SSF002", "SS019", "key13", "key200", "key108", "key28", "SS004", "key155", "ADF101", "SlogP_VSA9", "SSF015")
	
	testing_data = read.csv("Predictive_Model/model_v2/results_v2.4/topleft_15_1.3_2/independent_testing/BSEP_input.csv", header=TRUE, stringsAsFactors=F, sep=",", quote="", check.names=F)
	#feature = rbind(feature_data[c(pos_index,neg_index),6:ncol(feature_data)],testing_data[,4:ncol(testing_data)])
	
	testing_normalized = c()
	label_assign = c()
	label_predict = c()
	label_probability = c()
	for(i in 1:nrow(testing_data))
	{
	    feature = rbind(feature_data[c(pos_index,neg_index),selected_feature],testing_data[i,selected_feature])
	    feature = data.Normalization(feature,type="n5")  ### based all columns of feature to normalize
	    colnames(feature) = selected_feature
	    ## get training data: all pos and neg drugs and their features
	    ## define the columns to be selected
	    feature_training = feature[1:(length(pos_index)+length(neg_index)),selected_feature]
	    feature_training[is.na(feature_training)] = 0   ## replace na with 0, the reason not to remove na is to keep features consistant in training and testing
	    label_training = c( rep("pos", len=length(pos_index)), rep("neg", len=length(neg_index)) )
	    ##training_data = cbind(feature_data[c(pos_index,neg_index),c("drug_modeling", "fre_modeling", "logp_modeling", "RR_modeling", "year_modeling")], feature_training)
	    ##write.table(training_data, "Predictive_Model/model_v2/results_v2.4/independent_testing_2/training_data.txt", row.names=F, col.names=T, sep="\t", quote=F)
	    
		## get testing data
	    feature_testing = feature[(length(pos_index)+length(neg_index)+1):nrow(feature),selected_feature]  ## only extract selected feature
	    feature_testing[is.na(feature_testing)] = 0   ## replace na with 0, the reason not to remove na is to keep features consistant in training and testing
	    label_testing = testing_data[i,"label"]
	    drug_testing = testing_data[i,"ID"]
		temp = cbind(drug_testing, label_testing, feature_testing)
	    testing_normalized = rbind(testing_normalized, temp)
	    
		## predict
	    svm.model = svm(feature_training[,selected_feature], factor(label_training), cost=best_cost, kernel="linear", scale=F, cross=10, probability=TRUE)
	    svm.pred = predict(svm.model, feature_testing[,selected_feature], decision.values = T, probability=T) 
	
	    ## get true label and predict label, calculate performance
	    label_assign = c(label_assign, label_testing)  ## get the true label
	    temp_probability = as.vector(attr(svm.pred, "probabilities")[,"pos"])
		label_probability = c(label_probability, temp_probability) ## get the predicting probability of being pos
	    label_predict = c(label_predict, ifelse(temp_probability>=best_probability, "pos", "neg"))
	    
	}
	write.table(testing_normalized, "Predictive_Model/model_v2/results_v2.4/topleft_15_1.3_2/independent_testing/BSEP_data.txt", row.names=F, col.names=T, sep="\t", quote=F)
	mlabel = cbind(testing_data[,"ID"], label_assign, label_predict, label_probability)
	colnames(mlabel) = c("Drug", "Real", "Predict", "Probability")
	write.table(mlabel, "Predictive_Model/model_v2/results_v2.4/topleft_15_1.3_2/independent_testing/BSEP_label.txt", row.names=F, col.names=T, sep="\t", quote=F)
	
	### write the roc curve: only in the situation of knowing true label
	pdf(file="Predictive_Model/model_v2/results_v2.4/topleft_15_1.3_2/independent_testing/testing_roc.pdf", width=10, height=10)
	roc_curve = roc(label_assign, label_probability, levels=c("neg", "pos"))  ## levels help to define control(neg in this case) and case (pos in this case) samples in the data   
	## auc.polygon: whether fill in the roc area; grid: setup grid in the plot with interval; col: the color of roc curve
	plot(roc_curve, print.thres=T, print.auc=T, col="blue", auc.polygon=T, grid=c(0.1,0.2), grid.col=c("green", "red"), auc.polygon.col="lightcyan")
	t1 = paste("cost=", best_cost, sep="")
	t2 = paste("top=", best_top, sep="")
	text_value = c(t1, t2)
	legend("topleft", legend=text_value, text.font=4)  ## text.font: an integer specifying the font style of the legend text; possible values are: normal(1), bold(2), italic(3), bold and italic(4)
    dev.off()  
	

    ## write model summary, feature weights
	##write.table(capture.output(summary(svm.model)), "Predictive_Model/model_v2/results_v2.4/independent_testing_2/training_model_svm_summary.txt", row.names=F,  col.names=F, quote=F)  ## write the best model 
	##weight <- t(svm.model$coefs) %*% svm.model$SV   ## weight is a matrix: row=1, col=number of features
	##write.table(cbind(selected_feature, as.vector(weight)), "Predictive_Model/model_v2/results_v2.4/independent_testing_2/training_model_feature_weights.txt", row.names=T,  col.names=F, quote=F, sep="\t")
    
}

setwd("C:/Users/Desktop/")
##generate_drug_type ()
##preprocess()
##threshold_fre = 20   ## set up threshold for fre: assign labels to drugs indicating whether it is positive or negative based on their total number of reports
times <<- 500
for(t in 1:times)
{
    training_testing(t)
}
summarize_auc()
summarize_performance()
summarize_feature()



















