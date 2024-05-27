#The positive ion mode uses the following code

#Prepare the package to be used
library(xlsx)


File <- list.files(pattern=".xlsx")

if (length(grep("_Result",File))!=0){
  File <- File[-grep("_Result", File)]
}

#Read metabolic pathways, prototypes, metabolites tables
Route <- read.xlsx(File[5],sheetIndex = 1)
Route <- unique(Route)

Pro_drug <- read.xlsx(File[2],sheetIndex = 1)
Pro_drug <- subset(Pro_drug,!is.na(Pro_drug$RT))
Pro_drug <- subset(Pro_drug,!is.na(Pro_drug$HRMS..m.z))

Met_drug <- read.xlsx(File[4],sheetIndex = 1)
Met_drug <- subset(Met_drug,!is.na(Met_drug$RT))
Met_drug <- subset(Met_drug,!is.na(Met_drug$HRMS..m.z))

#Screen prototypes-------------------------------------------------------------------

Pro_drug_sample <- data.frame()#Extra prototypes after administration
Iso_drug_sample <- data.frame()#Extra isomer of prototypes after administration
Met_drug_screen <- data.frame()#Screen out prototypes

for (i in 1:length(Met_drug[,1])){
  
  #Determine whether or not it is prototype
  Judge1 <- subset(Pro_drug, 
                   Pro_drug$RT <= Met_drug$RT[i]+0.25 &
                     Pro_drug$RT >= Met_drug$RT[i]-0.25 &
                     Pro_drug$HRMS..m.z <= Met_drug$HRMS..m.z[i]*(1+5e-6) &
                     Pro_drug$HRMS..m.z >= Met_drug$HRMS..m.z[i]*(1-5e-6))
  
  if (nrow(Judge1) != 0) {
    
    New <- Met_drug[i,]
    New$Pro <- paste(Judge1$NO., collapse = " ")
    Pro_drug_sample <- rbind(Pro_drug_sample, New) 
    
  } else {
    Met_drug_screen <- rbind(Met_drug_screen, Met_drug[i,])
    
    Judge2 <- subset(Pro_drug, 
                     (Pro_drug$RT < Met_drug$RT[i]-0.5 |
                        Pro_drug$RT > Met_drug$RT[i]+0.5) &
                       Pro_drug$HRMS..m.z <= Met_drug$HRMS..m.z[i]*(1+5e-6) &
                       Pro_drug$HRMS..m.z >= Met_drug$HRMS..m.z[i]*(1-5e-6))
    
    if (nrow(Judge2) != 0) {
      
      New <- Met_drug[i,]
      New$Pro <- paste(Judge2$NO., collapse = " ")
      Iso_drug_sample <- rbind(Iso_drug_sample, New) 
      
    }
    
  }
  
  #Determining whether or not it is isomer of prototype
  
}

#Summarize the extra prototypes after administration and output to sheet1
write.xlsx(Pro_drug_sample, file = paste(File[4],"_Result.xlsx"),
           sheetName = "Prototypes",append=F)

write.xlsx(Iso_drug_sample, file = paste(File[4],"_Result.xlsx"),
           sheetName = "Isomers of prototype",append=T)

Met_drug_screen <- subset(Met_drug_screen,!Met_drug_screen$NO. %in% Iso_drug_sample$NO.)

#Screen one-route predictable metabolites of prototypes-------------------------------------------------

#Establishment of one-route predictable metabolites of prototypes list
Compound <- Pro_drug$NO.
One_Route_name <- Route[,1]
Add <- c("","[P+H]＆[M+NH4]","[P+NH4]＆[M+H]",
         "[P+H]＆[M+Na]","[P+Na]＆[M+H]",
         "[P+NH4]＆[M+Na]","[P+Na]＆[M+NH4]")

One_Route <- as.data.frame(outer(Route[,2], Pro_drug$HRMS..m.z, "+")) #Rows are metabolic pathways, columns are prototypes

One_Route_add <- list()
One_Route_add[[1]] <- One_Route
One_Route_add[[2]] <- One_Route + 14.003 + 1.008*3 #+(NH4-H)
One_Route_add[[3]] <- One_Route - 14.003 - 1.008*3 #-(NH4-H)
One_Route_add[[4]] <- One_Route + 22.990 - 1.008 #+(Na-H)
One_Route_add[[5]] <- One_Route - 22.990 + 1.008 #-(Na-H)
One_Route_add[[6]] <- One_Route + 22.990 - 14.003 - 1.008*4 #+(Na-NH4)
One_Route_add[[7]] <- One_Route - 22.990 + 14.003 + 1.008*4 #-(Na-NH4)

#Find and score 70 for one-route predictable metabolites
One_Route_pro <- data.frame() 

for (k in 1:length(Met_drug_screen[,1])){
  
  for(l in 1: length(One_Route_add)){
    
    indices <- which(One_Route_add[[l]] <= Met_drug_screen$HRMS..m.z[k]*(1+5e-6) 
                     & One_Route_add[[l]] >= Met_drug_screen$HRMS..m.z[k]*(1-5e-6),
                     arr.ind = TRUE) #5ppm
    indices <- as.data.frame(indices)
    
    if (length(indices[,1])!=0){
      
      for(m in 1:length(indices[,1])){
        
        Record <- Met_drug_screen[k,]  
        Record$Pro <- Compound[indices$col[m]]
        Record$Route <- One_Route_name[indices$row[m]]
        Record$Add <- Add[l]
        Record$FirMark <- 70 #Score 70
        
        One_Route_pro <- rbind(One_Route_pro, Record)
        
      }
    }
  }
}

#One-route predictable metabolites rely on MS2 to screen out-----------------------------------------

One_Route_pro$Dec <- One_Route_pro$HRMS..m.z - floor(One_Route_pro$HRMS..m.z)#Make sure it's 4 decimal places

#Circulating metabolites per signal
for (j in 1:length(One_Route_pro[,1])){
  
  Sec <- One_Route_pro$MS2[j]
  Sec1 <- strsplit(Sec,split = ", ")
  Sec1 <- unlist(Sec1)
  
  Int1 <- sub("(\\d+\\.).*", "\\1", Sec1)#Get integer part and decimal point,chr
  Dec1 <- sub("(\\d+)\\.(\\d+)\\(.*", "\\2", Sec1)#Get decimal part,chr
  Dec1 <- as.numeric(Dec1)
  
  #Score 0 if there is no match,score 20 if there is one match, and score 30 if there are two or more matches
  Match <- Pro_drug[Pro_drug$NO. == One_Route_pro$Pro[j],]
  
  if(nrow(Match) != 0){
    
    Sec <- Match$MS2
    Sec2 <- strsplit(Sec,split = ",")#Two files are split in different ways
    Sec2 <- unlist(Sec2)
    Int2 <- sub("(\\d+\\.).*", "\\1", Sec2)#Get integer part and decimal point,chr
    
    Same <- intersect(Int1,Int2)
    
    if(length(Same)==0){
      One_Route_pro$SecMark[j] <- 0
    }else if(length(Same) == 1) {
      One_Route_pro$SecMark[j] <- 20
    }else One_Route_pro$SecMark[j] <- 30
    
    #MS2 mass defects within ±25mDa, score 2 for each one
    addition <- sum(abs(Dec1-10000*One_Route_pro$Dec[j]) <= 500)
    One_Route_pro$MARK1[j] <- 2*addition #Note that the decimal portion requires 4
    
  } 
  
}

One_Route_pro <- subset(One_Route_pro, One_Route_pro$SecMark > 0)
One_Route_pro$SecMark <- One_Route_pro$SecMark + One_Route_pro$MARK1
One_Route_pro <- One_Route_pro[, -which(names(One_Route_pro) 
                                        %in% c("Dec","MARK1"))]

One_Route_pro$MARK <- One_Route_pro$FirMark + One_Route_pro$SecMark 
write.xlsx(One_Route_pro, file = paste(File[4],"_Result.xlsx"),
           sheetName = "One-route pre met of pro",append=T)

#Screening two-route predictable metabolites of prototypes-------------------------------------------------

#Establishment of two-route predictable metabolites of prototypes list
Compound <- Pro_drug$NO.
Route_name <- Route[,1]

Two_Route_ms <- as.data.frame(outer(Route[,2], Route[,2], "+"))
Two_Route_ms[lower.tri(Two_Route_ms)] <- -100000
Two_Route_ms <- unlist(Two_Route_ms)
Two_Route_name <- vector()

for (n in 1: length(Two_Route_ms)){
  Two_Route_name[n] <- paste(Route_name[ceiling(n/length(Route_name))],
                             Route_name[n - floor(n/length(Route_name)-0.0000000001)*length(Route_name)],
                             sep="＆") }

Two_Route <- data.frame(Two_Route_name,Two_Route_ms)
Two_Route <- subset(Two_Route, Two_Route$Two_Route_ms > -100000)#Completion of the list
Two_Route_name <- Two_Route[,1]

Add <- c("","[P+H]＆[M+NH4]","[P+NH4]＆[M+H]",
         "[P+H]＆[M+Na]","[P+Na]＆[M+H]",
         "[P+NH4]＆[M+Na]","[P+Na]＆[M+NH4]")

Two_Route_met <- as.data.frame(outer(Two_Route[,2], Pro_drug$HRMS..m.z, "+")) #Rows are metabolic pathways, columns are prototypes

Two_Route_add <- list()
Two_Route_add[[1]] <- Two_Route_met
Two_Route_add[[2]] <- Two_Route_met + 14.003 + 1.008*3 #+(NH4-H)
Two_Route_add[[3]] <- Two_Route_met - 14.003 - 1.008*3 #-(NH4-H)
Two_Route_add[[4]] <- Two_Route_met + 22.990 - 1.008 #+(Na-H)
Two_Route_add[[5]] <- Two_Route_met - 22.990 + 1.008 #-(Na-H)
Two_Route_add[[6]] <- Two_Route_met + 22.990 - 14.003 - 1.008*4 #+(Na-NH4)
Two_Route_add[[7]] <- Two_Route_met - 22.990 + 14.003 + 1.008*4 #-(Na-NH4)


#Find and score 70 for two-route predictable metabolites
Two_Route_pro <- data.frame() 

Met_drug_screen1 <- subset(Met_drug_screen, !Met_drug_screen$NO. %in% One_Route_pro$NO.)
#Met_drug_screen1 <- Met_drug_screen

for (k in 1:length(Met_drug_screen1[,1])){
  
  for(l in 1: length(Two_Route_add)){
    
    indices <- which(Two_Route_add[[l]] <= Met_drug_screen1$HRMS..m.z[k]*(1+5e-6) & 
                       Two_Route_add[[l]] >= Met_drug_screen1$HRMS..m.z[k]*(1-5e-6),
                     arr.ind = TRUE)
    indices <- as.data.frame(indices)
    
    if (length(indices[,1])!=0){
      
      for(m in 1:length(indices[,1])){
        
        Record <- Met_drug_screen1[k,]  
        Record$Pro <- Compound[indices$col[m]]
        Record$Route <- Two_Route_name[indices$row[m]]
        Record$Add <- Add[l]
        Record$FirMark <- 70 #Score 70
        
        Two_Route_pro <- rbind(Two_Route_pro, Record)
        
      }
    }
  }
}

#Two-route predictable metabolites rely on MS2 to screen out-----------------------------------------

Two_Route_pro$Dec <- Two_Route_pro$HRMS..m.z - floor(Two_Route_pro$HRMS..m.z)#Note that the decimal portion requires 4

#Circulating metabolites per signal
for (j in 1:length(Two_Route_pro[,1])){
  
  Sec <- Two_Route_pro$MS2[j]
  Sec1 <- strsplit(Sec,split = ", ")
  Sec1 <- unlist(Sec1)
  
  Int1 <- sub("(\\d+\\.).*", "\\1", Sec1)#Get integer part and decimal point,chr
  Dec1 <- sub("(\\d+)\\.(\\d+)\\(.*", "\\2", Sec1)#Get decimal part,chr
  Dec1 <- as.numeric(Dec1)
  
  #Score 0 if there is no match,score 20 if there is one match, and score 30 if there are two or more matches
  Match <- Pro_drug[Pro_drug$NO. == Two_Route_pro$Pro[j],]
  
  if(nrow(Match) != 0){
    
    Sec <- Match$MS2
    Sec2 <- strsplit(Sec,split = ",")#Two files are split in different ways
    Sec2 <- unlist(Sec2)
    Int2 <- sub("(\\d+\\.).*", "\\1", Sec2)#Get integer part and decimal point,chr
    
    Same <- intersect(Int1,Int2)
    
    if(length(Same)==0){
      Two_Route_pro$SecMark[j] <- 0
    }else if(length(Same) == 1) {
      Two_Route_pro$SecMark[j] <- 20
    }else Two_Route_pro$SecMark[j] <- 30
    
    #MS2 mass defects within ±25mDa, score 2 for each one
    addition <- sum(abs(Dec1-10000*Two_Route_pro$Dec[j]) <= 500)
    Two_Route_pro$MARK1[j] <- 2*addition #Note that the decimal portion requires 4
    
  } 
  
}

Two_Route_pro <- subset(Two_Route_pro, Two_Route_pro$SecMark > 0)
Two_Route_pro$SecMark <- Two_Route_pro$SecMark + Two_Route_pro$MARK1
Two_Route_pro <- Two_Route_pro[, -which(names(Two_Route_pro) 
                                        %in% c("Dec","MARK1"))]

Two_Route_pro$MARK <- Two_Route_pro$FirMark + Two_Route_pro$SecMark 
write.xlsx(Two_Route_pro, file = paste(File[4],"_Result.xlsx"),
           sheetName = "Two-route pre met of pro",append=T)

#Screening unpredictable metabolites of prototypes--------------------------------------

#Screen out predictable metabolites
Met_drug_screen_1 <- subset(Met_drug_screen, 
                            !Met_drug_screen$NO. %in% unique(One_Route_pro$NO.) &
                              !Met_drug_screen$NO. %in% unique(Two_Route_pro$NO.) )


#Initiate selection of unpredictable metabolites
if(length(Met_drug_screen_1[,1]) != 0){
  
  Met_drug_screen_1$MZInt <- as.integer(floor(Met_drug_screen_1$HRMS..m.z))
  Met_drug_screen_1$MZDec <- 
    Met_drug_screen_1$HRMS..m.z - floor(Met_drug_screen_1$HRMS..m.z)
  
  #Read the number of pre-screening data N elements
  Met_drug_screen_1$Nnum <- Met_drug_screen_1$Formula
  Met_drug_screen_1$Nnum <- gsub(".*(N\\d+).*","\\1",Met_drug_screen_1$Nnum)#Require a molecular formula beginning with C
  Met_drug_screen_1$Nnum <- gsub(".*( N ).*","\\1",Met_drug_screen_1$Nnum)
  Met_drug_screen_1$Nnum <- gsub(".*( N)$","\\1",Met_drug_screen_1$Nnum)
  Met_drug_screen_1$Nnum <- gsub(" ","",Met_drug_screen_1$Nnum)
  Met_drug_screen_1$Nnum <- gsub("C.*","0",Met_drug_screen_1$Nnum)
  Met_drug_screen_1$Nnum <- gsub("N(\\d+)","\\1",Met_drug_screen_1$Nnum)
  Met_drug_screen_1$Nnum <- gsub("N","1",Met_drug_screen_1$Nnum)
  Met_drug_screen_1$Nnum <- as.numeric(Met_drug_screen_1$Nnum)
  
  Pro_drug_1st <- Pro_drug
  Pro_drug_1st$MZInt <- as.integer(floor(Pro_drug_1st$HRMS..m.z))
  Pro_drug_1st$MZDec <- Pro_drug_1st$HRMS..m.z - floor(Pro_drug_1st$HRMS..m.z)
  
  #Read the number of pre-screening data N elements
  Pro_drug_1st$Nnum <- Pro_drug_1st$Formula
  Pro_drug_1st$Nnum <- gsub(".*(N\\d+).*","\\1",Pro_drug_1st$Nnum)#Require a molecular formula beginning with C
  Pro_drug_1st$Nnum <- gsub(".*( N ).*","\\1",Pro_drug_1st$Nnum)
  Pro_drug_1st$Nnum <- gsub(".*( N)$","\\1",Pro_drug_1st$Nnum)
  Pro_drug_1st$Nnum <- gsub(" ","",Pro_drug_1st$Nnum)
  Pro_drug_1st$Nnum <- gsub("C.*","0",Pro_drug_1st$Nnum)
  Pro_drug_1st$Nnum <- gsub("N(\\d+)","\\1",Pro_drug_1st$Nnum)
  Pro_drug_1st$Nnum <- gsub("N","1",Pro_drug_1st$Nnum)
  Pro_drug_1st$Nnum <- as.numeric(Pro_drug_1st$Nnum)
  
  Final_table1 <- data.frame()
  
  for(o in 1:length(Met_drug_screen_1[,1])){
    
    Judge3 <- subset(Pro_drug_1st, 
                     abs(Pro_drug_1st$MZInt - Met_drug_screen_1$MZInt[o]) < 100)
    Judge3 <- subset(Judge3, 
                     abs(Judge3$MZDec - Met_drug_screen_1$MZDec[o]) < 0.025)
    Judge3 <- subset(Judge3,
                     abs(Judge3$Nnum - Met_drug_screen_1$Nnum[o]) <= 3)
    
    if(nrow(Judge3) != 0){
      
      Addition <- Met_drug_screen_1[o,
                                    -((length(Met_drug_screen_1[1,])-2):length(Met_drug_screen_1[1,]))]
      
      Route <- paste("May be unpredictable metabolite of ",Judge3$NO.,sep = "")
      Route <- as.data.frame(Route)
      
      Addition <- tidyr::crossing(Addition, Route)
      Addition$FirMark <- 50
      
      Final_table1 <- rbind(Final_table1, Addition)
      
    }
  }
  
  #Unpredictable metabolites rely on MS2 to screen out-------------------------------------------
  
  if(nrow(Final_table1)!=0){
    
    Met_drug_2nd  <- subset(Final_table1, Final_table1$MS2 != "-" & 
                              !is.na(Final_table1$MS2))
    
    Met_drug_2nd$Pro <- gsub("May be unpredictable metabolite of (.*)", "\\1", Met_drug_2nd$Route)
    Met_drug_2nd$Dec2 <- Met_drug_2nd$HRMS..m.z - floor(Met_drug_2nd$HRMS..m.z)#Note that the decimal portion requires 4
    
    Pro_drug_2nd  <- subset(Pro_drug, Pro_drug$MS2 != "-" & !is.na(Pro_drug$MS2))
    
    #Circulating metabolites per signal
    for (j in 1:length(unlist(Met_drug_2nd[,1]))){
      
      Sec <- Met_drug_2nd$MS2[j]
      Sec1 <- strsplit(Sec,split = ", ")
      Sec1 <- unlist(Sec1)
      
      Int <- sub("(\\d+\\.).*", "\\1", Sec1)#Get integer part and decimal point,chr
      Dec <- sub("(\\d+)\\.(\\d+)\\(.*", "\\2", Sec1)#Get decimal part,chr
      Dec <- as.numeric(Dec)
      
      #Score 0 if there is no match,score 20 if there is one match, and score 30 if there are two or more matches
      Match <- Pro_drug_2nd[Pro_drug_2nd$NO. == Met_drug_2nd$Pro[j],]
      
      if(nrow(Match) != 0){
        
        Sec <- Match$MS2
        Sec2 <- strsplit(Sec,split = ",")#Two files are split in different ways
        Sec2 <- unlist(Sec2)
        Int2 <- sub("(\\d+\\.).*", "\\1", Sec2)#Get integer part and decimal point,chr
        
        Same <- intersect(Int,Int2)
        
        if(length(Same)==0){
          Met_drug_2nd$SecMark[j] <- 0
        }else if(length(Same) == 1) {
          Met_drug_2nd$SecMark[j] <- 20
        }else Met_drug_2nd$SecMark[j] <- 30
        
        #MS2 mass defects within ±25mDa, score 2 for each one
        addition <- sum(abs(Dec-10000*Met_drug_2nd$Dec2[j]) <= 250)
        Met_drug_2nd$MARK1[j] <- 2*addition #Note that the decimal portion requires 4
        
      } else{
        Met_drug_2nd_No <- rbind(Met_drug_2nd_No,Match)}
      
    }
    
    Met_drug_2nd <- subset(Met_drug_2nd, Met_drug_2nd$SecMark > 0)
    Met_drug_2nd$SecMark <- Met_drug_2nd$SecMark + Met_drug_2nd$MARK1
    Met_drug_2nd <- Met_drug_2nd[, -which(names(Met_drug_2nd) %in% c("Pro","Dec2","MARK1"))]
    Final_table1 <- Met_drug_2nd
    
  }
  
}

Final_table1$MARK <- Final_table1$FirMark + Final_table1$SecMark
Final_table1 <- subset(Final_table1, Final_table1$MARK > 50)

Met_drug_screen_2 <- subset(Met_drug_screen, 
                            !Met_drug_screen$NO. %in% unique(One_Route_pro$NO.) &
                              !Met_drug_screen$NO. %in% unique(Two_Route_pro$NO.) & 
                              !Met_drug_screen$NO. %in% unique(Final_table1$NO.)
)

if(nrow(Met_drug_screen_2) != 0){
  
  Met_drug_screen_2$FirMark <- 0
  Met_drug_screen_2$SecMark <- 0
  
  for(p in 1:nrow(Met_drug_screen_2)){
    Met_MS1 <- Met_drug_screen_2[p,3]
    Met_MS1_Int <- floor(Met_MS1)
    Met_MS1_Dec <- round(Met_MS1 - Met_MS1_Int, 4)*10^4
    
    Met_MS2 <- Met_drug_screen_2[p,6]
    Met_MS2 <- strsplit(Met_MS2, split = ", ")
    Met_MS2 <- unlist(Met_MS2)
    Met_MS2_Int <- sub("(\\d+)\\..*", "\\1", Met_MS2)
    Met_MS2_Int <- as.numeric(Met_MS2_Int)
    Met_MS2_Dec <- sub("\\d+\\.(\\d+).*", "\\1", Met_MS2)
    Met_MS2_Dec <- as.numeric(Met_MS2_Dec)
    Met_MS2 <- data.frame(Met_MS2_Int, Met_MS2_Dec)
    
    addition1 <- sum(abs(Met_MS2_Dec - Met_MS1_Dec) < 250)
    Met_drug_screen_2$SecMark[p] = 2*addition1
  }
  
  no_route_MDF_score <- data.frame()

  for(q in 1:nrow(Met_drug_screen_2)){

  Met_MS2 <- Met_drug_screen_2[q,6]
  Met_MS2 <- strsplit(Met_MS2, split = ", ")
  Met_MS2 <- unlist(Met_MS2)
  Met_MS2_Int <- sub("(\\d+)\\..*", "\\1", Met_MS2)
  Met_MS2_Int <- as.numeric(Met_MS2_Int)
  Met_MS2_Dec <- sub("\\d+\\.(\\d+).*", "\\1", Met_MS2)
  Met_MS2_Dec <- as.numeric(Met_MS2_Dec)
  Met_MS2 <- data.frame(Met_MS2_Int, Met_MS2_Dec)
  
  for(r in 1:nrow(Pro_drug)){
    Pro_name <- Pro_drug$NO.[r]
    Pro_MS2 <- Pro_drug[r,]$MS2
    Pro_MS2 <- unlist(strsplit(Pro_MS2, split = ","))
    Pro_MS2_Int <- sub("(\\d+)\\..*", "\\1", Pro_MS2)
    Pro_MS2_Int <- as.numeric(Pro_MS2_Int)
    Pro_MS2_Dec <- sub("\\d+\\.(\\d+).*", "\\1", Pro_MS2)
    if(nchar(Pro_MS2_Dec)[1] == 5){
      Pro_MS2_Dec = round(as.numeric(Pro_MS2_Dec)/10,0)
    }
    Pro_MS2_Dec <- as.numeric(Pro_MS2_Dec)
    Pro_MS2 <- data.frame(Pro_MS2_Int, Pro_MS2_Dec)
    
    get_MS2_MDF_count <- function(signal){
      data <- subset(Pro_MS2,abs(as.numeric(Pro_MS2_Int) - as.numeric(signal[1])) < 50)
      if(nrow(data)==0) return(0)
      count <- sum(abs(as.numeric(data$Pro_MS2_Dec) - as.numeric(signal[2])) <250, na.rm = T)
      return(count)
    }
    
    score <- sum(apply(Met_MS2, 1, get_MS2_MDF_count))*5
    no_route_MDF_score <- rbind(no_route_MDF_score,c(Met_drug_screen_2$NO.[q],Pro_name,score + Met_drug_screen_2$SecMark[q]))
    
  }
  
}
names(no_route_MDF_score) <- c("Met_Peak","Pro_Peak","MARK")

Final_table3 <- merge(Met_drug_screen_2,no_route_MDF_score, by.x = "NO.", by.y = "Met_Peak")
Final_table3 <- subset(Final_table3, as.numeric(Final_table3$MARK) > 50)
if(nrow(Final_table3 != 0)){
Final_table3$SecMark <- Final_table3$MARK 
Final_table3 <- merge(Final_table3,Met_drug_screen_1[,c("NO.","Nnum")],by = "NO.",all.x=T)
Final_table3 <- merge(Final_table3,Pro_drug_1st[,c("NO.","Nnum")],by.x="Pro_Peak", by.y = "NO.",all.x=T)
Final_table3 <- subset(Final_table3, abs(Final_table3$Nnum.y- Final_table3$Nnum.x)<=3)

Final_table3$Route <- paste("May be unpredictable metabolite of ",Final_table3$Pro_Peak,sep = "")

Final_table3 <- Final_table3[,names(Final_table1)]
Final_table1 <- rbind(Final_table1,Final_table3)
Final_table1 <- Final_table1[order(Final_table1$NO.),]
}


write.xlsx(Final_table1, file = paste(File[4],"_Result.xlsx"),
           sheetName = "Unpre met of pro", append=T)

}


if (length(Iso_drug_sample) != 0){
  
  Met_drug_screen_2 <- subset(Met_drug_screen, 
                            !Met_drug_screen$NO. %in% unique(One_Route_pro$NO.) &
                              !Met_drug_screen$NO. %in% unique(Two_Route_pro$NO.) & 
                              !Met_drug_screen$NO. %in% unique(Final_table1$NO.)
)

  Iso_drug_2nd  <- subset(Iso_drug_sample, 
                          Iso_drug_sample$MS2 != "-" & !is.na(Iso_drug_sample$MS2))
  
  Iso_screen <- data.frame()
  
  for (j in 1:length(Met_drug_screen_2[,1])){
    
    Sec <- Met_drug_screen_2$MS2[j]
    Sec1 <- strsplit(Sec,split = ", ")
    Sec1 <- unlist(Sec1)
    
    Int <- sub("(\\d+\\.).*", "\\1", Sec1)#Get integer part and decimal point,chr
    
    Same <- sapply(Int, 
                   function(element) grep(element, Iso_drug_2nd$MS2))
    Same <- unlist(Same)
    
    if (length(Same) != 0) {
      PASS <- Met_drug_screen_2[j,]
      PASS$Isomer <- paste(Iso_drug_2nd$NO.[unlist(Same)], collapse = " ")
      Iso_screen <- rbind(Iso_screen, PASS)
    }
    
  }

  write.xlsx(Iso_screen, file = paste(File[4],"_Result.xlsx"),
           sheetName = "Isomer_related_peak",append=T)                             
  
}

#No metabolic pathway identified---------------------------------------------------------------

#Score 0
Met_drug_screen_2 <- subset(Met_drug_screen, 
                            !Met_drug_screen$NO. %in% unique(One_Route_pro$NO.) &
                              !Met_drug_screen$NO. %in% unique(Two_Route_pro$NO.) & 
                              !Met_drug_screen$NO. %in% unique(Final_table1$NO.) &
                              !Met_drug_screen$NO. %in% unique(Iso_screen$NO.)
)

Met_drug_screen_2 <- subset(Met_drug_screen_2, 
                              !Met_drug_screen_2$NO. %in% Iso_screen$NO.)

if(nrow(Met_drug_screen_2) != 0){
  
  Met_drug_screen_2$Route <- "No route match"
  Met_drug_screen_2$FirMark <- 0
  Met_drug_screen_2$SecMark <- 0
  
  for(p in 1:nrow(Met_drug_screen_2)){
    Met_MS1 <- Met_drug_screen_2[p,3]
    Met_MS1_Int <- floor(Met_MS1)
    Met_MS1_Dec <- round(Met_MS1 - Met_MS1_Int, 4)*10^4
    
    Met_MS2 <- Met_drug_screen_2[p,6]
    Met_MS2 <- strsplit(Met_MS2, split = ", ")
    Met_MS2 <- unlist(Met_MS2)
    Met_MS2_Int <- sub("(\\d+)\\..*", "\\1", Met_MS2)
    Met_MS2_Int <- as.numeric(Met_MS2_Int)
    Met_MS2_Dec <- sub("\\d+\\.(\\d+).*", "\\1", Met_MS2)
    Met_MS2_Dec <- as.numeric(Met_MS2_Dec)
    Met_MS2 <- data.frame(Met_MS2_Int, Met_MS2_Dec)
    
    addition1 <- sum(abs(Met_MS2_Dec - Met_MS1_Dec) < 250)
    Met_drug_screen_2$SecMark[p] = 2*addition1
  }
    
  Met_drug_screen_2$MARK <- Met_drug_screen_2$FirMark + Met_drug_screen_2$SecMark
  write.xlsx(Met_drug_screen_2, file = paste(File[4],"_Result.xlsx"),
             sheetName = "No route", append=T)
  
}