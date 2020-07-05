#in order, a= g salt 1, b= mL acid to titrate salt 1, c= g salt 2, d= mL acid for salt 2
fullcalculation<-function(a,b,c,d){

library(dplyr)
library(readxl)
library(tidyr)

#periodic table data
ptable<-read_excel("Data/elements.xls")

#this will be replaced with inputs for the Shiny app
sampledata<-tibble(sample_number=c(1,2), mass_g=c(a,c), titration_volume_mL=c(b, d))

#the next two formulas simply create strings for general formulas based on metal charge
carbonateformulapicker<-function(x){
if(x==1){
  formula<-"M2CO3"
} else if (x ==2){
  formula<-"MCO3"
} else if (x ==3){
  formula<-"M2(CO3)3"
} else if (x==4){
  formula<-"M(CO3)2"
} else{
  formula<-"none"
}
  formula
}

hydroxideformulapicker<-function(x){
  if(x==1){
    formula<-"MOH"
  } else if (x ==2){
    formula<-"M(OH)2"
  } else if (x ==3){
    formula<-"M(OH)3"
  } else if (x==4){
    formula<-"M(OH)4"
  } else{
    formula<-"none"
  }
  formula
}

#the next two formulas calculate the acid molarity that would lead to a given molar mass based on sample data
carbonatemolaritycalc<-function(molarmass,valence,volume,samplemass){
  if(valence ==1){
    molarity<-1/((molarmass*2+60.008)*(10*volume/(2*1000*samplemass)))
  } else if (valence ==2){
    molarity<-1/((molarmass+60.008)*(10*volume/(2*1000*samplemass)))
  } else if (valence==3){
    molarity<-1/((molarmass+60.008*2)*(10*volume/(4*1000*samplemass)))
  } else if (valence==4){
    molarity<-1/((molarmass*2+60.008*3)*(10*volume/(6*1000*samplemass)))
  } else{
    molarity<-"none"
  }
  molarity
}

hydroxidemolaritycalc<-function(molarmass,valence,volume,samplemass){
  if(valence ==1){
    molarity<-1/((molarmass+17.069)*(10*volume/(1000*samplemass)))
  } else if (valence ==2){
    molarity<-1/((molarmass+17.069*2)*(10*volume/(2*1000*samplemass)))
  } else if (valence==3){
    molarity<-1/((molarmass+17.069*3)*(10*volume/(3*1000*samplemass)))
  } else if (valence==4){
    molarity<-1/((molarmass+17.069*4)*(10*volume/(4*1000*samplemass)))
  } else{
    molarity<-"none"
  }
  molarity
}

#this is a wrapper function that selects one of the two above formulas depending on the identity of the salt
molaritypicker<-function(molarmass,valence,salt_type,volume,samplemass){
  if (grepl("CO3",salt_type,fixed=TRUE)){
    carbonatemolaritycalc(molarmass,valence,volume,samplemass)
  }
  else if(grepl("OH",salt_type,fixed=TRUE)){
    hydroxidemolaritycalc(molarmass,valence,volume,samplemass)
  }
  else{
    molarity<-"none"
    molarity
  }
}

#the following two formulas are the reverse of the pair above, it calculates a molar mass given a molarity and sample data
carbonatemmcalc<-function(molarity,valence,volume,samplemass){
  if(valence ==1){
    molarmass<-(samplemass/10/(molarity*(volume/1000)/2)-60.008)/2
  } else if (valence ==2){
    molarmass<-(samplemass/10/(molarity*(volume/1000)/2)-60.008)
  } else if (valence==3){
    molarmass<-(samplemass/10/(molarity*(volume/1000)/6)-(60.008*3))/2
  } else if (valence==4){
    molarmass<-(samplemass/10/(molarity*(volume/1000)/4)-(60.008*2))
  } else{
    molarmass<-"none"
  }
  molarmass
}

hydroxidemmcalc<-function(molarity,valence,volume,samplemass){
  if(valence ==1){
    molarmass<-(samplemass/10/(molarity*(volume/1000))-17.0069)
  } else if (valence ==2){
    molarmass<-(samplemass/10/(molarity*(volume/1000)/2)-(17.0069*2))
  } else if (valence==3){
    molarmass<-(samplemass/10/(molarity*(volume/1000)/3)-(17.0069*3))
  } else if (valence==4){
    molarmass<-(samplemass/10/(molarity*(volume/1000)/4)-(17.0069*4))
  } else{
    molarmass<-"none"
  }
  molarmass
}

#this wrapper function picks one of the two formulas above based on salt identity
mmpicker<-function(valence,salt_type, volume,samplemass, molarity){
  if (grepl("CO3",salt_type,fixed = TRUE)){
    carbonatemmcalc(molarity,valence,volume,samplemass)
  }
  else if(grepl("OH",salt_type,fixed = TRUE)){
    hydroxidemmcalc(molarity,valence,volume,samplemass)
  }
  else{
    molarmass<-"none"
    molarmass
  }
}
#taking only what is needed from the periodic table and creates one line per charge per metal
metalstable<-ptable%>%filter(!Group %in% c("Non-Metal","Halogen","Noble Gas"))%>%select(Sym,"Atomic Weight",Valence)%>%separate_rows( Valence, sep=",")%>%transform(Valence=as.numeric(Valence))%>%transform(Atomic.Weight=abs(Atomic.Weight))%>%filter(Valence<5, Valence>0)

#every metal+charge combination appears twice here, once as a carbonate, once as a hydroxide, salts_list is the basis for all further calc
carbonate_formulas<-metalstable%>%mutate(saltFormula=sapply(metalstable$Valence,carbonateformulapicker))
hydroxide_formulas<-metalstable%>%mutate(saltFormula=sapply(metalstable$Valence,hydroxideformulapicker))
salts_list<-merge(carbonate_formulas,hydroxide_formulas, all=TRUE)

#call to molaritypicker to create a version of salts_list with virtual molarities attached
salt1molarities<-salts_list%>%mutate(salt_1_molarity= as.numeric(mapply(molaritypicker,valence=salts_list$Valence, molarmass=salts_list$Atomic.Weight, salt_type=salts_list$saltFormula, MoreArgs=list(sampledata[1,3],sampledata[1,2]))))

#this takes a given molarity and finds the salt with the smallest difference between the predicted and actual metal molar mass
salt2calc<- function(molarity){
  salt2table<-salts_list%>%mutate(calcmm=as.numeric(mapply(mmpicker, valence=salts_list$Valence, salt_type=salts_list$saltFormula, MoreArgs=list(sampledata[2,3],sampledata[2,2], molarity))),mmdiff=abs(calcmm-salts_list$Atomic.Weight))%>%arrange(mmdiff)
  salt2table[1,]%>%mutate(salt_1_molarity=molarity)
}

#this takes every molarity calculated above and finds the best paired metal based on sample data, since the test is done with two groups at a time, this takes a long time to run
salt2<-as.data.frame(t(sapply(salt1molarities$salt_1_molarity,salt2calc)))

#this combines the salt 1 and 2 results, finds a virtual molarity halfway between the ideal one for each salt, then returns results in order from the lowest molar mass error
#there is also some output cleanup here
#this will be the product for the Shiny app
finalresults<-merge(salt1molarities,salt2, by="salt_1_molarity", all=FALSE)%>%mutate(salt_2_molarity= as.numeric(mapply(molaritypicker,valence=Valence.y, molarmass=Atomic.Weight.y, salt_type=saltFormula.y, MoreArgs=list(sampledata[2,3],sampledata[2,2]))),Molarity=(salt_1_molarity+salt_2_molarity)/2,CalculatedMolarMassMetal1=mapply(mmpicker, valence=Valence.x, salt_type=saltFormula.x, molarity=Molarity, MoreArgs=list(sampledata[1,3],sampledata[1,2])),mmdiffavg=(abs(as.numeric(CalculatedMolarMassMetal1)-as.numeric(Atomic.Weight.x)))+as.numeric(mmdiff)/2, CalculatedMolarMassMetal2=mapply(mmpicker, valence=Valence.y, salt_type=saltFormula.y, molarity=Molarity, MoreArgs=list(sampledata[2,3],sampledata[2,2])))%>%arrange(as.numeric(mmdiffavg))%>%select(Molarity,Sym.x,Atomic.Weight.x,saltFormula.x,CalculatedMolarMassMetal1,Sym.y,Atomic.Weight.y,saltFormula.y,CalculatedMolarMassMetal2)%>%rename(Metal1=Sym.x,MolarMassMetal1=Atomic.Weight.x,Formula1=saltFormula.x,Metal2=Sym.y,MolarMassMetal2=Atomic.Weight.y,Formula2=saltFormula.y)
finalresults[1:40,]
}
