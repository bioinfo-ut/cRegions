# Created by Mikk Puustusmaa
# Last edit on 25.07.2017

#------------------------------SETTINGS-------------------------------------
#mode = 1 # Single positions
#mode = 2 # sliding window mode
#mode = 3 # both (web tool uses mode 3)
args <- commandArgs(TRUE)
mode = as.numeric(args[1])


#A metric is not calculated for a column in MSA which has more gaps than allowed
allowedGap = as.double(as.double(args[2])/100)

#In sliding window mode a column in MSA is skipped during arithmetic mean calculation instead of terminating if there is more gaps than skipGap allows.
#TODO set limit for skipping
skipGap = as.double(as.double(args[3])/100)


#-----sliding window mode-----

#window size - How many values in arithmetic mean calculation
windowSize = as.numeric(args[4])

# Position in three-nucleotide codon which is used in arithmetic mean calculation. 
# E.g. if codonPosition = 3 and window size = 3, then arithmetic mean is calculated over 3,6,9 etc.
codonPosition = as.numeric(args[5])

#location of data
folder = args[6]
observedDataFile = toString(args[7]) #count data
observedProportionDataFile = toString(args[8]) #observed nucleotide proportions 
predictedUniformFile = toString(args[9]) #weighted or normal
predictedFile = toString(args[10]) #weighted or normal
codonUsageFile = toString(args[11])


#zeroProbabilityCorrection = FALSE
#if(args[12] == 1){zeroProbabilityCorrection = TRUE}
#zeroProbabilityCorrection = TRUE

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


# codonUsage.txt -> info: 	49	2073
numberSeq  = as.numeric(sapply(strsplit(firstLine, "\t",fixed = FALSE),"[[", 2))
aligmentSoftware = sapply(strsplit(firstLine, "\t",fixed = FALSE),"[[", 3)

#Total length of MSA including gaps
N = length(observed$pos)


colorPos1 = "blue"
colorPos2 = "green"
colorPos3 = "red"

#Array for large table
rStatistics=seq(1, N, by=1)

#Array for large table
rStatistics_sliding_window=seq(1, N, by=1)


createPlotChisq <- function(values, title){
  fileName = paste(folder,"/",title,".txt", sep="")
  location = paste(folder,"/",title, sep="")
  
  write.table(values,file = fileName, row.names = FALSE, col.names = FALSE, sep="\t")
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
  fileConn<-file(paste(location,"_avg.txt", sep=""))
  writeLines(jsArray, fileConn)
  close(fileConn)
  
  values <- sapply(values, as.numeric)
  rStatistics_sliding_window <<- cbind(rStatistics_sliding_window,values) # also in createPlot
}

createPlot <- function(values,pos1,pos2,pos3, title, currentMode){

	if(currentMode == 1 || currentMode == 3){
		#File location
		fileName = paste(folder,"/",title,".txt", sep="")
		
		#Write table to file
		write.table(values,file = fileName, row.names = FALSE, col.names = FALSE, sep="\t")	

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
		
		for (i in 1:N ) {
		  if (title=="CHISQ" || title=="UNIFORM_CHISQ") {
				
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


		location = paste(folder,"/",title, sep="")
		

		jsArray1 = paste(jsArray1, "]", sep="")
		fileConn<-file(paste(location,"_pos1.txt",sep=""))
		writeLines(jsArray1, fileConn)
		close(fileConn)

		jsArray2 = paste(jsArray2, "]", sep="")
		fileConn<-file(paste(location,"_pos2.txt",sep=""))
		writeLines(jsArray2, fileConn)
		close(fileConn)

		jsArray3 = paste(jsArray3, "]", sep="")
		fileConn<-file(paste(location,"_pos3.txt",sep=""))
		writeLines(jsArray3, fileConn)
		close(fileConn)

		#-----------------------
		#For highcharts
		#----------------------

	}

if(currentMode == 2 || currentMode == 3){
		addedValues = rep(NA,N)
		yMax = 0


		#-----------------------
		#For highcharts
		#----------------------
		
		jsArray = "["
		for(i in seq(codonPosition, (length(values)-3*windowSize), by = 3)){
			sum = 0
			#If in some positions in MSA, there is insertion that is caused only by one or few sequences, 
			#it is not good to terminate arithmetic mean calculation, rather skip that position
			counter_PositionsInSum = 0
			currentPosition = i
			

			#If there is a position (column) in MSA, where is more gaps than allowed,
			#then this positioni is set to NA and if there is NA in sliding window
			#value for this window is not calculated
			NAinWindow = FALSE
			
		
			val = values[currentPosition]
			repeat{
				if (is.na(val)){
					if (counter_PositionsInSum == 0){
						#first position
						NAinWindow = TRUE
						break
					}
					else if(observed[currentPosition,6] < skipGap){
						NAinWindow = TRUE
						break
						
					}
				}else{
				  val = as.double(val)
					sum = sum + val
					counter_PositionsInSum = counter_PositionsInSum + 1
					
				}
				
				currentPosition = currentPosition + 3

				
  				if(counter_PositionsInSum == windowSize){
    					break
  				}else if (currentPosition >= length(values)-3){
						NAinWindow = TRUE
						break
					
				}


			}#Repeat END
	
			if(!NAinWindow){
				addedValues[i] = (sum / windowSize)
				value = round((sum / windowSize), digits = 3)
				jsArray = paste(jsArray, paste("[",i,",",value,"], ", sep=""), sep="")
				
				if (yMax < (sum / windowSize)){
					yMax = (sum / windowSize)
				}
			}#NA if end
		

	}#for end

		if (title!="CHISQ" & title!="UNIFORM_CHISQ") {
			xMax = length(addedValues)

			location = paste(folder,"/",title, sep="")

			jsArray = paste(jsArray, "]", sep="")
			fileConn<-file(paste(location,"_mean.txt", sep=""))
			writeLines(jsArray, fileConn)
			close(fileConn)
			
			#-----------------------
			#Mean
			#----------------------
			mean = mean(addedValues, na.rm=TRUE) 
			addedValuesSD = sd(addedValues, na.rm=TRUE)
			
			writeToFile = c(mean, addedValuesSD)
			location = paste(folder,"/",title, sep="")
			fileConn = file(paste(location,"_arithmeticMean_ofAverage.txt", sep=""))
				write(writeToFile, fileConn)
			close(fileConn)
		}
		

		#-----------------------
		#For highcharts
		#----------------------

		addedValues <- sapply(addedValues, as.numeric)
		#data into one table
		if (is.null(rStatistics_sliding_window )){
			if (title!="CHISQ" & title!="UNIFORM_CHISQ") {
				rStatistics_sliding_window<<-addedValues
			} 
		}
		else{
			
			if(title!="CHISQ" & title!="UNIFORM_CHISQ"){
				rStatistics_sliding_window <<- cbind(rStatistics_sliding_window,addedValues) # also in createPlotChisq
			}
		}
	
}
  
	
}


calculateRMSD <- function(observed,predicted){

	difA = (observed[1] - predicted[1])^2
	difC = (observed[2] - predicted[2])^2
	difG = (observed[3] - predicted[3])^2
	difT = (observed[4] - predicted[4])^2
	val = sqrt((difA + difC + difG + difT)/4)
	return (val);
	
}

calculateMAXDIF <- function(observed,predicted){

	difA = abs(observed[1] - predicted[1])
	difC = abs(observed[2] - predicted[2])
	difG = abs(observed[3] - predicted[3])
	difT = abs(observed[4] - predicted[4])
	v = c(difA, difC, difG, difT)
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
	values = rep(NA,N)
	pos1 = rep(NA,N)
	pos2 = rep(NA,N)
	pos3 = rep(NA,N)
	
	for(i in 1:N) {
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
	
	createPlot(values,pos1,pos2,pos3, metric,mode)
}



#------------------------------REAL Pearson's Chi-squared Test for Count Data for sliding window mode------------------------------
#


#pos	 A 	 C 	 G 	 T 	ctrl	gap
#10	0.3333	0.5385	0.0769	0.0513	1.0	0.2041
chisqTestValues = rep(NA,N)

for(i in seq(codonPosition, (length(chisqTestValues)-3*windowSize), by = 3)){
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

createPlotChisq(chisqTestValues,"CHISQ")

colnames(rStatistics) <- c("MSA_pos","RMSD","RMSD(uniform codon usage)","MAXDIF","MAXDIF(uniform codon usage)","CHISQ p-value","CHISQ p-value (uniform codon usage)")
dataName = paste(folder,"/raw_values.tsv", sep="")
write.table(rStatistics,file = dataName, row.names = FALSE, col.names = TRUE, sep="\t")


colnames(rStatistics_sliding_window) <- c("MSA_pos","RMSD","RMSD(uniform codon usage)","MAXDIF","MAXDIF(uniform codon usage)","CHISQ p-value")
dataName = paste(folder,"/raw_values_sliding_window.tsv", sep="")
write.table(rStatistics_sliding_window,file = dataName, row.names = FALSE, col.names = TRUE, sep="\t")


