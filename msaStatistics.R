# Created by Mikk Puustusmaa
# Last edit on 17.04.2018

#------------------------------SETTINGS-------------------------------------
args <- commandArgs(TRUE)

#A metric is not calculated for a column in MSA which has more gaps than allowed
allowedGap = as.double(as.double(args[1])/100)

#In sliding window mode a column in MSA is skipped during arithmetic mean calculation instead of terminating if there is more gaps than skipGap allows.
#TODO set limit for skipping
skipGap = as.double(as.double(args[2])/100)


#-----sliding window mode-----

#window size - How many values are included in the sliding window
windowSize = as.numeric(args[3])

# Position in three-nucleotide codon which is used in the arithmetic mean calculation. 
# E.g. if codonPosition = 3 and window size = 3, then arithmetic mean is calculated over 3,6,9 etc.
codonPosition = as.numeric(args[4])

#location of data
folder = args[5]
observedDataFile = toString(args[6]) #count data
observedProportionDataFile = toString(args[7]) #observed nucleotide proportions 
predictedUniformFile = toString(args[8]) #weighted or normal
predictedFile = toString(args[9]) #weighted or normal
codonUsageFile = toString(args[10])


#------------------------------SETTINGS END--------------------------------
location = paste(folder, "/", observedDataFile, sep="")
observedCount = read.table(location,sep="\t",dec=".",header=T)


location = paste(folder, "/", observedProportionDataFile, sep="")
observed = read.table(location,sep="\t",dec=".",header=T)


location = paste(folder, "/", predictedUniformFile, sep="")
predictedUniform = read.table(location,sep="\t",dec=".",header=T)


location = paste(folder, "/", predictedFile, sep="")
predicted  = read.table(location,sep="\t",dec=".",header=T)

location = paste(folder, "/", codonUsageFile, sep="")
con = file(location, "rt") 

firstLine = readLines(con, 1) # Read one line
arr <- strsplit(firstLine, '\t') [[1]]
numberSeqArr <- strsplit(arr[1], ':') [[1]]


# codon_usage_bias.tsv -> number_of_seq: 49	 msa_length:2073
numberSeq <- as.numeric(numberSeqArr[2])

#Total length of MSA including gaps
codonAlignmentLength = length(observed$pos)


colorPos1 = "blue"
colorPos2 = "green"
colorPos3 = "red"

rStatistics=seq(1, codonAlignmentLength, by=1)
rStatistics_sliding_window=seq(1, codonAlignmentLength, by=1)


createChartDataChisq <- function(values, title){
  resultFile = paste(folder,"/",title,".txt", sep="")
  location = paste(folder,"/",title, sep="")
  
  write.table(values,file = resultFile, row.names = FALSE, col.names = FALSE, sep="\t")
  jsArray = "["
  for (i in 1:length(values) ) {
    pvalue = values[i]
	if(!is.na(pvalue) && pvalue == 0){pvalue = 1e-300}
	val = (-1)*log10(pvalue)
	if(!is.na(val) && val == "Inf"){val = 1e-300}
	
    if (!is.na(val)  && !is.null(val)){
      val = as.double(val)
      val = round(val, digits = 3)
      jsArray = paste(jsArray, paste("[",i,",",val,"], ", sep=""), sep="")
    }
  }
  
  jsArray = paste(jsArray, "]", sep="")
  fileConn<-file(paste(location,"_sliding_window_chart.txt", sep=""))
  writeLines(jsArray, fileConn)
  close(fileConn)
  
  values <- sapply(values, as.numeric)
  rStatistics_sliding_window <<- cbind(rStatistics_sliding_window,values) 
}

createSinglePositionChartData <- function(singlePositionValues,pos1,pos2,pos3,metric){
  #File location
  resultFile = paste(folder,"/",metric,".txt", sep="")
  
  #Write table to file
  write.table(values,file = resultFile, row.names = FALSE, col.names = FALSE, sep="\t")	

  values <- sapply(values, as.numeric)
  if (is.null(rStatistics)){
    rStatistics<<-values
  }
  else{
    rStatistics <<- cbind(rStatistics,values)
  }
  #-----------------------
  #For highcharts
  #----------------------
  jsArray1 = "["
  jsArray2 = "["
  jsArray3 = "["
  
  for (i in 1:codonAlignmentLength ) {
    if (metric=="CHISQ" || metric=="UNIFORM_CHISQ") {
      
      pvalue = values[i]
      if(!is.na(pvalue) && pvalue == 0){pvalue = 1e-50}
      val = (-1)*log10(pvalue)
      if(!is.na(val) && val == "Inf"){val = 1e-300}
    }else {
      val = values[i]
    }
    
    if (!is.na(val) && !is.null(val)){
      val = as.double(val)
      val = round(val, digits = 3)
      if( i%%3 ==1 ){
        jsArray1 = paste(jsArray1, paste("[",i,",", val,"], ", sep=""), sep="")
      }
      else if( i%%3 == 2 ){
        jsArray2 = paste(jsArray2, paste("[",i,",", val,"], ", sep=""), sep="")
      }
      else{
        jsArray3 = paste(jsArray3, paste("[",i,",", val,"], ", sep=""), sep="")
      }
    }
    
  }
  location = paste(folder,"/",metric, sep="")
  jsArray1 = paste(jsArray1, "]", sep="")
  fileConn<-file(paste(location,"_pos1_chart.txt",sep=""))
  writeLines(jsArray1, fileConn)
  close(fileConn)

  jsArray2 = paste(jsArray2, "]", sep="")
  fileConn<-file(paste(location,"_pos2_chart.txt",sep=""))
  writeLines(jsArray2, fileConn)
  close(fileConn)

  jsArray3 = paste(jsArray3, "]", sep="")
  fileConn<-file(paste(location,"_pos3_chart.txt",sep=""))
  writeLines(jsArray3, fileConn)
  close(fileConn)
}

createSlidingWindowChartData <- function(singlePositionValues,windowSize,codonPosition,metric){
  slidingWindowValues = rep(NA,codonAlignmentLength)
  jsArray = "["
  for(i in seq(codonPosition, (length(values)-3*windowSize), by = 3)){
    sum = 0
    counter_PositionsInSum = 0
    currentPosition = i
    NAinWindow = FALSE
   
    for(k in seq(currentPosition, (length(values)-3*windowSize), by = 3)){
      val = singlePositionValues[k]
      gapRatio = observed[k,6]
      if (gapRatio >= skipGap){next}
      else if (is.na(val)){NAinWindow = TRUE;break}
      else {
        val = as.double(val)
        sum = sum + val
        counter_PositionsInSum = counter_PositionsInSum + 1
        currentPosition = currentPosition + 3
      }

      if (counter_PositionsInSum == windowSize){
          break
      } else if (currentPosition >= length(values)-3){
        NAinWindow = TRUE
        break
      }
    }#Repeat END

    if(!NAinWindow){
      avg = sum / windowSize
      slidingWindowValues[i] = (avg)
      value = round(avg, digits = 3)
      jsArray = paste(jsArray, paste("[",i,",",value,"], ", sep=""), sep="")
    }#NA if end
  

}#for end

  if (metric!="CHISQ" & metric!="UNIFORM_CHISQ") {
    #Write sliding window values to file
    location = paste(folder,"/",metric, sep="")
    jsArray = paste(jsArray, "]", sep="")
    fileConn<-file(paste(location,"_sliding_window_chart.txt", sep=""))
    writeLines(jsArray, fileConn)
    close(fileConn)
  
    mean = mean(slidingWindowValues, na.rm=TRUE) 
    slidingWindowValuesSD = sd(slidingWindowValues, na.rm=TRUE)
    writeToFile = c(mean, slidingWindowValuesSD)
    location = paste(folder,"/",metric, sep="")
    fileConn = file(paste(location,"_sliding_window_mean_and_standard_deviation.txt", sep=""))
      write(writeToFile, fileConn)
    close(fileConn)
  }
  

  #-----------------------
  #For highcharts
  #----------------------

  slidingWindowValues <- sapply(slidingWindowValues, as.numeric)
  #data into one table
  if (is.null(rStatistics_sliding_window )){
    if (metric!="CHISQ" & metric!="UNIFORM_CHISQ") {
      rStatistics_sliding_window<<-slidingWindowValues
    } 
  }
  else{
    
    if(metric!="CHISQ" & metric!="UNIFORM_CHISQ"){
      rStatistics_sliding_window <<- cbind(rStatistics_sliding_window,slidingWindowValues) # also in createChartDataChisq
    }
  }
	
}


calculateRMSD <- function(observed,predicted){

	val = NA
	dif = 0
	predictedLength <- length(predicted)
	for(i in 1:predictedLength) {
		dif = dif + (observed[i] - predicted[i])^2
	}
	val = sqrt(dif/predictedLength)
	return (val);
	
}

calculateMAXDIF <- function(observed,predicted){

	val = NA
	v = c()
	predictedLength <- length(predicted)
	for(i in 1:predictedLength) {
		dif = abs(observed[i] - predicted[i])
		v = c(v,dif)
	}

	val = max(v)
	return (val);
	
}

calculateCHISQ <- function(observed,predicted){

	if(length(predicted)>1){
		test = chisq.test(observed,p = predicted, rescale.p=TRUE)
		pvalue = test$p.value
		val = pvalue
	} else {
		val = NA 
	}
	
	return (val);
	
}





metrics = c("RMSD","UNIFORM_RMSD","MAXDIF","UNIFORM_MAXDIF","CHISQ","UNIFORM_CHISQ")

for(metric in metrics){
	#NA is generally interpreted as a missing value and has various forms - NA_integer_, NA_real_, etc
	values = rep(NA,codonAlignmentLength)
	pos1 = rep(NA,codonAlignmentLength)
	pos2 = rep(NA,codonAlignmentLength)
	pos3 = rep(NA,codonAlignmentLength)
	
	for(i in 1:codonAlignmentLength) {
		if(observed[i,6] >= allowedGap){
			values[i] = NA
		}else{
			observedOnePosition = as.numeric(observed[i,2:5]) #observed values proportions
			observedCountOnePosition = as.numeric(observedCount[i,2:5])

			#-----------------BEGIN RMSD & MAXDIF-------------------------
			if (metric == "UNIFORM_RMSD" || metric == "UNIFORM_MAXDIF" || metric == "UNIFORM_CHISQ"){
				predictedOnePosition = as.numeric(predictedUniform[i,2:5]) # predicted values, uniform codon usage
			} else {
				predictedOnePosition = as.numeric(predicted[i,2:5]) # predicted values
			}
			
			#if a position only has an amino acid (e.g. tryptophan) with single codon, metric should not be calculated
			x = c()
			z = c()
			p = c() 
			for(j in 1:4) {
				if (predictedOnePosition[j] > 0){
					x[length(x)+1] = as.numeric(observedOnePosition[j])
					z[length(z)+1] = as.numeric(observedCountOnePosition[j])
					p[length(p)+1] = predictedOnePosition[j]
				} 
			}

			if (metric == "RMSD" || metric == "UNIFORM_RMSD") {
				val = calculateRMSD(x,p)
			} else if (metric == "MAXDIF" || metric == "UNIFORM_MAXDIF"){
				val = calculateMAXDIF(x,p)
			}

			
			#-----------------END RMSD & MAXDIF-------------------------
			
			#-----------------BEGIN  Pearson's Chi-squared Test for Count Data-------------------------
			#chisq.test performs chi-squared contingency table tests and goodness-of-fit tests
			
			if (metric == "CHISQ" || metric == "UNIFORM_CHISQ") {
				val = calculateCHISQ(z,p)
				
			}
			
			#-----------------END CHISQ-------------------------
			
			
			values[i] = val
			if(i%%3 == 0){pos3[i]= val}
			else if (i%%2 == 0){pos2[i]= val}
			else {pos1[i]= val}
		}
	}
	
	createSinglePositionChartData(values,pos1,pos2,pos3, metric)
  createSlidingWindowChartData(values,windowSize,codonPosition,metric)
}



#------------------------------REAL Pearson's Chi-squared Test for Count Data for sliding window mode------------------------------
#pos	 A 	 C 	 G 	 T 	ctrl	gap
#10	0.3333	0.5385	0.0769	0.0513	1.0	0.2041
chisqTestValues = rep(NA,codonAlignmentLength)

for(i in seq(codonPosition, (codonAlignmentLength-3*windowSize), by = 3)){
	if(observed[i,6] >= allowedGap){
		chisqTestValues[i] = NA
	}else{
	
		x = c()
		p = c()

		for(j in 2:5) {
		  if (predicted[i,j] != "NA" && observedCount[i,j] != "NA"){
				observed_count = as.numeric(observedCount[i,j])
				x[length(x)+1] = observed_count
				p[length(p)+1] = as.double(predicted[i,j])
		  }
		}
		
		count = 1
		isEnd = TRUE #If skipping reaches to the end of the MSA
		for(k in seq((i+3), length(chisqTestValues), by = 3)){
			if (count == windowSize) {isEnd = FALSE; break}
		
			gapRatio = observed[k,6]
			if(gapRatio >= skipGap){next}
		  
			if(gapRatio < allowedGap){
				for(j in 1:4) {
					x[j] = x[j] + as.numeric(observedCount[k,j+1])
					p[j] = p[j] + as.double(predicted[k,j+1])
				}
				count = count + 1
			}
		}


		if(isEnd){
			break # p-value is not calculated
		}else{
			window_x = round(x/windowSize)
			window_p = p/windowSize
			
			x = c()
			p = c()
			
			for(z in seq(4)) {
				if (window_p[z] > 0){
					x[length(x)+1] = window_x[z]
					p[length(p)+1] = window_p[z]
					#p[length(p)+1] = ifelse(window_p[z] > 0, as.double(window_p[z]), 1e-10)
				} 
			}

			if(length(x)>1){
			  test = chisq.test(x,p = p, rescale.p=TRUE)
			  pvalue = test$p.value
			  chisqTestValues[i] = pvalue

			}else{
			  chisqTestValues[i] = NA
			}
		}
	}
}

createChartDataChisq(chisqTestValues,"CHISQ")

colnames(rStatistics) <- c("MSA_pos","RMSD","RMSD(uniform codon usage)","MAXDIF","MAXDIF(uniform codon usage)","CHISQ p-value","CHISQ p-value (uniform codon usage)")
dataName = paste(folder,"/raw_values.tsv", sep="")
write.table(rStatistics,file = dataName, row.names = FALSE, col.names = TRUE, sep="\t")


colnames(rStatistics_sliding_window) <- c("MSA_pos","RMSD","RMSD(uniform codon usage)","MAXDIF","MAXDIF(uniform codon usage)","CHISQ p-value")
dataName = paste(folder,"/raw_values_sliding_window.tsv", sep="")
write.table(rStatistics_sliding_window,file = dataName, row.names = FALSE, col.names = TRUE, sep="\t")


