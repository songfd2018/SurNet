library(stringr) # String manipulation, Regex
library(ggplot2)
library(xtable)

# Parse the header of each email file
parseHeader <- function(header){
  #MessageID <- str_sub(str_extract(header, "Message-ID:.*"), start = 12)
  Date <- str_sub(str_extract(header,"Date:.*"), start = 7)
  
  #Date<-trans_time(Date)
  Date<-str_sub(Date,end=-7)
  #Conversion of dates
  ## UTC time
  datesTest <- strptime(Date, format = "%a, %d %b %Y %H:%M:%S %z", tz = "GMT")
  ## localtime
  datesLocal <- strptime(Date, format = "%a, %d %b %Y %H:%M:%S")
  
  From <- str_sub(str_extract(header,"X-From:.*"), start = 9)
  To <- str_sub(str_extract(header,"X-To:.*"), start = 7)
  Subject <- str_sub(str_extract(header,"Subject:.*"), start = 10)
  Cc <- str_sub(str_extract(header,"X-cc:.*"), start = 7)
  #X-bcc <- str_sub(str_extract(header,"X\\-bcc:.*"), start = 8)
  
  headerParsed <- data.frame(Date, datesTest, datesLocal, From, To, Cc, Subject,# X-bcc, 
                             stringsAsFactors = FALSE)
  return(headerParsed)
}

# Parse the body of each email file
parseBody <- function(body){
  num<-length(body)
  content<-rep(NA,num)
  forwarded<-rep(NA,num)
  id.forward<-rep(0,num)
  #forward.start<- str_locate(body, "---------------------- Forwarded by")[,1]#find the first forward benchmark
  #forward.start1<-str_locate(body, "----------------------")[,1]#find the first forward benchmark
  forward.start<- str_locate(body, "Forwarded by")[,1]#find the first forward benchmark
  for(i in 1:num){#excluding the forwarded information
    if(is.na(forward.start[i])){
      content[i]<-body[i]
    }else{
      content[i]<-str_sub(body[i],end = forward.start[i]-1)
      
      #collecting the forwarded information
      forwarded[i]<-str_sub(body[i],start = forward.start[i])
    }
  }
  id.forward[!is.na(forwarded)]<-1
  
  #deal with direct reply in the emails
  sub_in_content<-str_locate_all(content,"Subject")
  outlier<-NULL
  original<-rep(NA,num)
  content.end<-rep(NA,num)
  cat.bench<-rep(NA,num)
  time.start<-rep(NA,num)
  id.reply<-rep(0,num)
  original.date<-rep(NA,num)
  
  #content can be null
  for(i in 1:num){
    if(dim(sub_in_content[[i]])[1]>0){#The email is a directed reply email
      #\n\n\n\n\n and --------- are two possible ends of email
      poss.end<-str_locate(content[i],c("\n\n\n\n","\n\t\n\t\n","Original Message","you wrote:"))[,1]
      if(sum(!is.na(poss.end))>0){
        cat.bench[i]<-max(which(!is.na(poss.end)))
        content.end[i]<-poss.end[cat.bench[i]]
        if(cat.bench[i]<3){
          original[i]<-str_sub(content[i],start = content.end[i])
          content[i]<-str_sub(content[i],end=content.end[i]-1)
        }else if(cat.bench[i]==3){
          original[i]<-str_sub(content[i],start = content.end[i]-5)
          content[i]<-str_sub(content[i],end=content.end[i]-6)
          # #identify whether this emal is reply email or forwarding email
          # original.from<-str_sub(str_extract(original[i],"From:.*"), start = 7)
          # if(str_detect(original.from,"Crenshaw")|str_detect(original.from,"crenshaw")) id.reply[i]<-1
          # if(str_detect(original.from,"Gibne")|str_detect(original.from,"gibne")) id.reply[i]<-1
          # if(str_detect(original.from,"Kamins")|str_detect(original.from,"kamins")) id.reply[i]<-1
          # #extract sent, which is sending time of original email
          # if(id.reply==1){
          #   original.date[i]<-str_sub(str_extract(original[i],"Sent:.*"), start = 7)
          # }
        }else{
          original[i]<-str_sub(content[i],start = content.end[i]-29)
          content[i]<-str_sub(content[i],end=content.end[i]-30)
        }
        #extract time of original emails
        #time.start[i]<-str_match(original[i],"[0-9]")
        #time.locate<-str_locate(original[i],time.start)
        id.reply[!is.na(original)]<-1
        
      }else{
        outlier<-c(outlier,i)#record the outlier content to find the content
      }
    }
  }
  
  index.original<-which(!is.na(original))
  
  res<-data.frame(id.reply,#whether the mail has an original email
                  content,#email content
                  original,#original email
                  id.forward,#whether the mail forwards other email
                  forwarded,#fowarded message
                  stringsAsFactors = FALSE)
  return(res)
}

# parse the body of each file as an event
parseEvent<-function(header1,body){
  
  #output a list: each element is a chain of emails called an event, refered as a matrix
  #for each emails we need collect email ID time from to cc, content and relation to the preivous emails(reply, forward) 
  res<-NULL
  
  #regard each chain of emails as an event, 
  #first to count the number of emails in each event by counting the number of subjects
  loc.subject<-gregexpr("Subject:",body)
  num.event<-length(loc.subject)
  num.email<-rep(NA,num.event)
  for(i in 1:num.event){
    if(loc.subject[[i]][1]>0){
      num.email[i]<-length(loc.subject[[i]])+1
    }else{
      num.email[i]<-1
    }
  }
  
  #collect outlier index
  outlier<-NULL
  
    #identify each from, to, subject and content for each email
    for(i in 1:num.event){
      if(num.email[i]>1){
        temp1<-NULL #save the part of content before subject: body+header
        temp2<-body[i]#save the remaining part: subject+body
        temp2<-str_replace_all(temp2,"\n>>","\n")
        temp2<-str_replace_all(temp2,"\n>","\n")
        
        temp3<-matrix(NA,num.email[i],9)#event matrix
        colnames(temp3)<-c("Date","UCTtime","LOCALtime","From","To","Cc","Subject","State","Content")
        #temp3[,1:3] three date
        #temp3[,4] from
        #temp3[,5] to
        #temp3[,6] cc
        #temp3[,7] subject
        #temp3[,8] state: Single, End, Forwarded, Replied
        temp3[1,1:7]<-as.matrix(header1[i,])
        temp3[1,8]<-"end"
        repeat.subject<-NULL
        remaining<-NULL#to record the remaining part if there is false subject
        #loc.end<-diff(c(1,loc.subject[[i]])) #the location of subject 
        for(j in 2:num.email[i]){
          #find the first content
          if(is.na(temp2)){
            outlier<-c(outlier,i)
            break
          }
          loc.end<-regexpr("Subject:",temp2)
          temp1<-paste(remaining,str_sub(temp2,end = loc.end-1),sep="")
          temp2<-str_sub(temp2,start = loc.end)
          sub<-str_extract(temp2,"Subject:.*")
          Subject <- str_sub(sub, start = 9)
          temp2<-str_trim(str_sub(temp2,start=str_length(Subject)+9))
          
          if(loc.end>5){#Maybe repeated subject, if loc.end<5
          temp3[j,7]<-str_trim(Subject)
          
          #divide body and header
          benchmark.end<-c("Date:","To:","From:",
                           "@ENRON",".com>",
                           "\n\n\n\n","\n\t\n\t",
                           "Reply Separator",
                           "Forwarded by","Forwarded message",
                           "Original Message",
                           "Message-ID:","Message-Id")
          con.end<-str_locate(temp1,benchmark.end)[,1]#"you wrote:", "Forwarded by"))[,1]
          if(sum(!is.na(con.end))>0){
            cat.bench<-max(which(!is.na(con.end)))
            content.end<-con.end[cat.bench]
            if(cat.bench<6&min(con.end,na.rm=T)<content.end){
              content.end<-min(con.end,na.rm=T)
              cat.bench<-which(con.end==content.end)
            }
            if(cat.bench==1){
              header.back<-str_sub(temp1,start = content.end)
              temp3[j-1,9]<-str_sub(temp1,end=content.end-1)
              
              Date <- str_sub(str_extract(header.back,"Date:.*"),start=7)
              temp3[j,1]<-str_trim(Date)
              temp3[j,8]<-"OnlyDate"
            }#"Date"
            if(cat.bench==2){
              header.back<-str_sub(temp1,start = content.end)
              temp3[j-1,9]<-str_sub(temp1,end=content.end-1)
              
              To <- str_sub(str_extract(header.back,"To:.*"),start=5)
              Cc <- str_sub(str_extract(header.back,"cc:.*"),start=5)
              
              temp3[j,5:6]<-str_trim(c(To,Cc))
              temp3[j,8]<-"RepliedNoFromDate"
            }#"To:"
            if(cat.bench==3){
              header.back<-str_sub(temp1,start = content.end)
              temp3[j-1,9]<-str_sub(temp1,end=content.end-1)
              
              if(str_detect(header.back," on ")){
                From <- str_sub(str_extract(header.back,"From:.*"),start=7)
                Date <- strsplit(From," on ")[[1]][2]
                From <- strsplit(From," on ")[[1]][1]
              }
              if(str_detect(header.back,"Sent:")){
                From <- str_sub(str_extract(header.back,"From:.*"),start=7)
                Date <- str_sub(str_extract(header.back,"Sent:.*"),start=7)
              }
              
              To <- str_sub(str_extract(header.back,"To:.*"),start=5)
              Cc <- str_sub(str_extract(header.back,"Cc:.*"),start=5)
              
              temp3[j,1]<-str_trim(Date)
              temp3[j,4:6]<-str_trim(c(From,To,Cc))
              temp3[j,8]<-"Replied"
            }#"From:"
            if(cat.bench==4|cat.bench==5){
              body.rough<-str_sub(temp1,end = content.end)
              loc.bienter<-gregexpr("\n\n",body.rough)[[1]]
              content.end<-loc.bienter[length(loc.bienter)]
              if(content.end>=0){
                header.back<-str_sub(temp1,start = content.end)
                temp3[j-1,9]<-str_sub(temp1,end=content.end-1)
              }else{
                header.back<-temp1
                temp3[j-1,9]<-""
              }
              
              loc.date<-regexpr('\\d',header.back)
              loc.to<-regexpr('To:',header.back)
              loc.cc<-regexpr('cc:',header.back)
              From <- substr(header.back,1,loc.date-1)
              Date <- substr(header.back,loc.date,loc.date+19)
              To <- substr(header.back,loc.to+4,loc.cc-1)
              Cc <- substring(header.back,loc.cc+4)
              temp3[j,1]<-str_trim(Date)
              temp3[j,4:6]<-str_trim(c(From,To,Cc))
              temp3[j,8]<-"Replied"
              
            }#"@ENRON",".com>",
            if(cat.bench==6|cat.bench==7){
              header.back<-str_sub(temp1,start = content.end)
              temp3[j-1,9]<-str_sub(temp1,end=content.end-1)
              
              if(str_detect(header.back,".com> on")){
                loc.date<-regexpr('\\d',str_sub(temp1,start = con.end[5]))+con.end[5]-con.end[cat.bench]
              }else{
                loc.date<-regexpr('\\d',header.back)
              }
              loc.to<-regexpr('To:',header.back)
              loc.cc<-regexpr('cc:',header.back)
              From <- substr(header.back,1,loc.date-1)
              Date <- str_sub(header.back,start=loc.date)
              if(str_detect(Date,"AM")|str_detect(Date,"PM")){
                Date <-str_sub(Date,end=str_locate(Date,"M")[1])
              }else{
                Date <-str_split(Date,"\n",n=2,simplify = T)[1]
              }
              To <- substr(header.back,loc.to+4,loc.cc-1)
              Cc <- substring(header.back,loc.cc+4)
              temp3[j,1]<-str_trim(Date)
              temp3[j,4:6]<-str_trim(c(From,To,Cc))
              temp3[j,8]<-"Replied"
            
              }#"\n\n\n\n","\n\t\n\t",
            if(cat.bench==8){
              header.back<-str_split(temp2,"\n\n\n",n=2,simplify = T)[1]
              temp2<-str_split(temp2,"\n\n\n",n=2,simplify = T)[2]
              
              temp3[j-1,9]<-str_sub(temp1,end=content.end-1)
              
              From <- str_sub(str_extract(header.back,"Author:.*"),start=9)
              Date <- str_sub(str_extract(header.back,"Date:.*"),start=7)
              
              temp3[j,1]<-str_trim(Date)
              temp3[j,4]<-str_trim(From)
              temp3[j,8]<-"RepliedNoFromDate"
            }#"Reply Separator"
            if(cat.bench==9){
              header.back<-str_sub(temp1,start = content.end)
              temp3[j-1,9]<-str_sub(temp1,end=content.end-1)
              
              From<-NA
              Date<-NA
              
              if(str_detect(header.back,"\n\n\n")){
                header.back<-str_split(header.back,"\n\n\n",n=2,simplify = T)[2]
                loc.date<-regexpr('\\d',header.back)
                From <- substr(header.back,1,loc.date-4)
                
                Date <- str_sub(header.back,start=loc.date)
                Date <-str_split(Date,"\n",n=2,simplify = T)[1]
              }
              if(str_detect(header.back,"\n\t\n\t")){
                header.back<-str_split(header.back,"\n\t\n\t",n=2,simplify = T)[2]
                From <- str_sub(str_extract(header.back,"From:.*"),start=7)
                Date <- str_split(From,"          ",n=2,simplify = T)[2]
                From <- str_split(From,"          ",n=2,simplify = T)[1]
              }
              
              loc.to<-regexpr('To:',header.back)
              loc.cc<-regexpr('cc:',header.back)
              To <- substr(header.back,loc.to+4,loc.cc-1)
              Cc <- substring(header.back,loc.cc+4)
              temp3[j,1]<-str_trim(Date)
              temp3[j,4:6]<-str_trim(c(From,To,Cc))
              temp3[j,8]<-"Forwarded"
          
              }#"Forwarded by"
            if(cat.bench==10){
                header.back<-str_sub(temp1,start = content.end)
                temp3[j-1,9]<-str_sub(temp1,end=content.end-1)
                
                From <- str_sub(str_extract(header.back,"From:.*"),start=7)
                To <- str_sub(str_extract(header.back,"To:.*"),start=5)
                Cc <- str_sub(str_extract(header.back,"Cc:.*"),start=5)
                Date <- str_sub(str_extract(header.back,"Date:.*"),start=7)
                temp3[j,1]<-str_trim(Date)
                temp3[j,4:6]<-str_trim(c(From,To,Cc))
                temp3[j,8]<-"Forwarded"
                
              }#"Forwarded message",
            if(cat.bench==11){
              header.back<-str_sub(temp1,start = content.end)
              temp3[j-1,9]<-str_sub(temp1,end=content.end-1)
              
              From <- str_sub(str_extract(header.back,"From:.*"),start=7)
              To <- str_sub(str_extract(header.back,"To:.*"),start=5)
              Cc <- str_sub(str_extract(header.back,"Cc:.*"),start=5)
              Date<-NA
              if(str_detect(header.back,"Date:")){
                Date <- str_sub(str_extract(header.back,"Date:.*"),start=7)
              }
              if(str_detect(header.back,"Sent:")){
                Date <- str_sub(str_extract(header.back,"Sent:.*"),start=7)
              }
              temp3[j,1]<-str_trim(Date)
              temp3[j,4:6]<-str_trim(c(From,To,Cc))
              temp3[j,8]<-"Original"
            
              }#"Original Message",
            if(cat.bench==12){
              loc.from<-str_locate(temp1,"From:")[1]
              header.back<-str_sub(temp1,start = loc.from)
              temp3[j-1,9]<-str_sub(temp1,end= loc.from-1)
              
              From <- str_sub(str_extract(header.back,"From:.*"),start=7)
              To <- str_sub(str_extract(header.back,"To:.*"),start=5)
              Cc <- str_sub(str_extract(header.back,"Cc:.*"),start=5)
              Date <- str_sub(str_extract(header.back,"Date:.*"),start=7)
              temp3[j,1]<-str_trim(Date)
              temp3[j,4:6]<-str_trim(c(From,To,Cc))
              temp3[j,8]<-"Return-Path"
            }#"Message-ID:",
            if(cat.bench==13){
              header.back<-str_sub(temp1,start = content.end)
              temp3[j-1,9]<-str_sub(temp1,end=content.end-1)
              
              From <- str_sub(str_extract(header.back,"From:.*"),start=7)
              To <- str_sub(str_extract(header.back,"To:.*"),start=5)
              Cc <- str_sub(str_extract(header.back,"Cc:.*"),start=5)
              Date <- str_sub(str_extract(header.back,"Date:.*"),start=7)
              temp3[j,1]<-str_trim(Date)
              temp3[j,4:6]<-str_trim(c(From,To,Cc))
              temp3[j,8]<-"Return-Path"
            }#"Message-Id"
            
            remaining<-NULL
            
          }else{#This subject is part of 
            repeat.subject<-c(repeat.subject,j)
            remaining<-paste(temp1,sub,sep="")
            
            #outlier<-rbind(outlier,c(i,j))#record j-th email in i-th is an outlier.
          }
          
          }else{
            repeat.subject<-c(repeat.subject,j)
            remaining<-paste(temp1,sub,sep="")
          }
          
          
        }#end of for loop for emails
        temp3[j,9]<-str_trim(temp2)
        if(!is.null(repeat.subject)){
          temp3<-temp3[-repeat.subject,]
        }
        res[[i]]<-temp3
  
      }else{
        res[[i]]<-c(as.matrix(header1[i,]),"Single",body[i])
      }
      #iteratively extract the email
      #con.end<-str_locate(body[i],c("\n\n\n\n","\n\t\n\t\n","Original Message","you wrote:","Forwarded by"))[,1]
    }#end of events loop
  
  return(res)
  
}

unique_email <- function(time,usr){
  
  args=list(Time = time, User = usr) #iteration setting
  
  res<-1
  index<-1
  temp.time <- time[1]
  temp.usrname <- usr[1]
  n<-length(time)
  for(i in 1:n){
    t<-time[i]
    usrname<-usr[i]
    if(t!=temp.time|usrname!=temp.usrname){
      res<-c(res,index)
      temp.time<-t
      temp.usrname<-usrname
    }
    index<-index+1
  }
  
  return(res)
}

index_no_repeat <- function(time,usr){
  
  res<-1
  index<-1
  temp.time <- time[1]
  temp.usrname <- usr[1]
  n<-length(time)
  for(i in 1:n){
    t<-time[i]
    usrname<-usr[i]
    if(t!=temp.time|usrname!=temp.usrname){
      res<-c(res,index)
      temp.time<-t
      temp.usrname<-usrname
    }
    index<-index+1
  }
  
  return(res)
}

#exclude no-letter and no-number chars in the string
prep_fun <- function(x) {
  x %>% 
    # make text lower case
    str_to_lower %>% 
    # remove non-alphanumeric symbols
    str_replace_all("[^[:alpha:]]", " ") %>% 
    # collapse multiple spaces
    str_replace_all("\\s+", " ")
}

#deal with the problem of not figuring out letters in the strptime function
trans_time <- function(time){
  #transform the Date to standard form
  time<-str_replace(time,"Monday","1")
  time<-str_replace(time,"Tuesday","2")
  time<-str_replace(time,"Wednesday","3")
  time<-str_replace(time,"Thursday","4")
  time<-str_replace(time,"Friday","5")
  time<-str_replace(time,"Saturday","6")
  time<-str_replace(time,"Sunday","0")
  
  time<-str_replace(time,"Mon","1")
  time<-str_replace(time,"Tue","2")
  time<-str_replace(time,"Wed","3")
  time<-str_replace(time,"Thu","4")
  time<-str_replace(time,"Fri","5")
  time<-str_replace(time,"Sat","6")
  time<-str_replace(time,"Sun","0")
  
  time<-str_replace(time,"January","01")
  time<-str_replace(time,"February","02")
  time<-str_replace(time,"March","03")
  time<-str_replace(time,"April","04")
  time<-str_replace(time,"May","05")
  time<-str_replace(time,"June","06")
  time<-str_replace(time,"July","07")
  time<-str_replace(time,"August","08")
  time<-str_replace(time,"September","09")
  time<-str_replace(time,"October","10")
  time<-str_replace(time,"November","11")
  time<-str_replace(time,"December","12")
  
  time<-str_replace(time,"Jan","01")
  time<-str_replace(time,"Feb","02")
  time<-str_replace(time,"Mar","03")
  time<-str_replace(time,"Apr","04")
  time<-str_replace(time,"May","05")
  time<-str_replace(time,"Jun","06")
  time<-str_replace(time,"Jul","07")
  time<-str_replace(time,"Aug","08")
  time<-str_replace(time,"Sep","09")
  time<-str_replace(time,"Oct","10")
  time<-str_replace(time,"Nov","11")
  time<-str_replace(time,"Dec","12")
  
  #Polish
  time<-str_replace(time,"lutego","02")
  
  
  return(time)
}

#unify the intra-event sending time
unify_time <- function(time){
  res<-NA
  rpt<-0#detect double time accepatance
  time<-trans_time(time)
  
  if(!is.na(strptime(time, "%w, %d %m %Y %H:%M"))){#Weekday, Day Month Year Hour:Min
    res<-strptime(time, "%w, %d %m %Y %H:%M")
    rpt<-rpt+1
  }
  
  if(!is.na(strptime(time,"%w, %m %d, %Y %H:%M"))){#Weekday, Month Day, Year hour(1-12):Min
    AP<-str_sub(time,start=-2)
    res<-strptime(time,"%w, %m %d, %Y %H:%M")
    if(AP=="PM"){
      res<-res+12*3600
    }
    rpt<-rpt+1
    #time.intra[j]<-reply.time-send.time
  }
  
  if(!is.na(strptime(time,"%d %m %Y %H:%M"))){#Day Month Year Hour:Min
    AP<-str_sub(time,start=-2)
    res<-strptime(time,"%d %m %Y %H:%M")
    if(AP=="PM"){
      res<-res+12*3600
    }
    rpt<-rpt+1
  }
  
  if(!is.na(strptime(time, "%d/%m/%Y %H:%M"))){#Day/Month/Year(1999) Hour:Min
    AP<-str_sub(time,start=-2)
    res<-strptime(time,"%d/%m/%Y %H:%M")
    if(AP=="PM"){
      res<-res+12*3600
    }
    rpt<-rpt+1
  }
  
  if(!is.na(strptime(time, "%m/%d/%y %H:%M"))){#Month/Day/year(99) Hour:Min
    AP<-str_sub(time,start=-2)
    res<-strptime(time,"%m/%d/%y %H:%M")
    if(res$year<99) res$year=res$year+1900
    if(res$year<99) res$year=res$year+100
    if(AP=="PM"){
      res<-res+12*3600
    }
    rpt<-rpt+1
  }
  
  if(!is.na(strptime(time, "%m/%d/%Y %H:%M"))){#Month/Day/Year(1999) Hour:Min
    AP<-str_sub(time,start=-2)
    res<-strptime(time,"%m/%d/%Y %H:%M")
    if(AP=="PM"){
      res<-res+12*3600
    }
    rpt<-rpt+1
  }
  
  if(!is.na(strptime(time, "%y-%m-%d %H:%M"))){#year(99)-Month-Day/ Hour:Min
    AP<-str_sub(time,start=-2)
    res<-strptime(time,"%y-%m-%d %H:%M")
    if(res$year<99) res$year=res$year+1900
    if(res$year<99) res$year=res$year+100
    if(AP=="PM"){
      res<-res+12*3600
    }
    rpt<-rpt+1
  }
  
  if(rpt>1){
    res<--1
  }
  return(res)
}

#data preprocessing for Enron email dataset
#read emails.csv
data<-read.csv("RawData/emails.csv", stringsAsFactors = FALSE)

#file split 
fileSplit <- str_split(data$file, "/")

###################################################
# Extract the folder names which emails belong to #
###################################################
n<-nrow(data)
leng.filename<-rep(NA,n)
user<-rep(NA,n)
for(i in 1:n){
  #leng.filename[i]<-length(fileSplit[[i]])
  user[i]<-fileSplit[[i]][1]
}
user.list<-table(user)

####################################
#extract the content from the email#
####################################

breaks <- str_locate(data$message, "\n\n")
headers <- str_sub(data$message, end = breaks[,1] - 1 )
bodies <- str_sub(data$message, start = breaks[,2] + 1)

###################
#Parse the headers#
###################

prev <- Sys.getlocale("LC_TIME"); Sys.setlocale("LC_TIME", "C")
headerParsed<- parseHeader(headers)

########################################
#sort by time and exclude repeat emails#
########################################

datesTest<-as.POSIXlt(headerParsed$datesTest)
#select time window 99-02
datesTest$year<-(datesTest$year+10)%%100+90
events.interest<-which(datesTest$year>98&datesTest$year<103)
headerParsed$datesTest<-datesTest
headerParsed<-headerParsed[events.interest,]
emails <- data$message[events.interest]
bodies <- bodies[events.interest]
user <- user[events.interest]

#order from the earliest to latest
order.time <- order(headerParsed$datesTest)
emails <- data$message[order.time]
headerParsed <- headerParsed[order.time,]
bodies <- bodies[order.time]
user <- user[order.time]

# #index of no-repeat
# index.nrep<-index_no_repeat(headerParsed.usr1$datesTest)
# index.nrep<-index.nrep[-1]#the first time is a outlier value 

index.nrep <- index_no_repeat(headerParsed$datesTest,user)

emails <- emails[index.nrep]
headerParsed <- headerParsed[index.nrep,]
bodies <- bodies[index.nrep]
user <- user[index.nrep]

#free variables
rm(data)
rm(breaks)
rm(fileSplit)

########################
#Parse chains of events#
########################
events<-parseEvent(headerParsed,bodies)

#################################################################################
# Extract the events whose sender and receiver both belong to the employee list #
#################################################################################
#extract content
n.event<-length(events)

##########################################
# Mutually select the users' identifiers #

#select sending emails
employees<-read.table("RawData/employeelist.txt",stringsAsFactors = FALSE, header = T)
n.usr<-nrow(employees)


# load the user identifers
user.iden <- read.csv(file = "RawData/user_identifier.csv", header = TRUE, stringsAsFactors = FALSE)
identifier_list <- list()
for(usr in 1:n.usr){
  temp <- unlist(user.iden[usr,3:6])
  identifier_list[[usr]] <- tolower(temp[which(temp!="")])
}


##########################################################################
# Extract events whose sender and receiver are both in the employee list #
##########################################################################
interest_record <- NULL
# To record: 1. Event index
#            2. Folder name
#            3. Sender
#            4. Receiver
#            5. Local Time
#            6. UTC Time
#            7. Subject
#            8. Category
#            9. Content
count <- 0
for(evt in 1:n.event){
  temp.event <- events[[evt]]
  n.emails <- nrow(temp.event)

  # for each email, search the sender and the receiver in the employee list, respectively
  if(!is.matrix(temp.event)){
    temp.event[4:6] <- tolower(temp.event[4:6])
    
    sender_no <- NULL
    receiver_no <- NULL
    for(usr in 1:n.usr){
      # get the identifier of the "usr"
      temp_id <- identifier_list[[usr]]
      if(all(!is.na(temp.event[4:6]))){
        if(sum(str_detect(temp.event[4],temp_id))>0){
          sender_no <- c(sender_no, usr)
        }
        if(sum(str_detect(temp.event[5],temp_id)) + sum(str_detect(temp.event[6],temp_id)) >0){
          receiver_no <- c(receiver_no, usr)
        }
      }
    }
    if(!is.null(sender_no) & !is.null(receiver_no) & !identical(sender_no,receiver_no)){
      # count the number of words in the email content
      content <- temp.event[9]
      content <- gsub("[^[:alnum:]]", " ", x= content)
      num_words <- length(unlist(strsplit(x=content, split=" +")))
      if(length(receiver_no) == 1){
        evt_record <- c(evt, # Event index,
                        user[evt], # folder name
                        sender_no, # sender
                        employees[sender_no, 4], # name of sender
                        receiver_no, # receiver
                        employees[receiver_no,4], # name of receiver
                        temp.event[1], # Original time
                        temp.event[2], # UTC time
                        temp.event[3], # Local time
                        length(receiver_no), # Numeber of receiver
                        gsub("\""," ",temp.event[7]), # subject
                        temp.event[8], # category
                        num_words, # number of words in the last email
                        num_words, # number of words in the original email
                        gsub("\n|\""," ",temp.event[9])) # content
      }else{
        evt_record <- cbind(evt, # Event index,
                            user[evt], # folder name
                            sender_no, # sender
                            employees[sender_no, 4], # name of sender
                            receiver_no, # receiver
                            employees[receiver_no,4], # name of receiver
                            temp.event[1], # Original time
                            temp.event[2], # UTC time
                            temp.event[3], # Local time
                            length(receiver_no), # Numeber of receiver
                            gsub("\""," ",temp.event[7]), # subject
                            temp.event[8], # category
                            num_words, # number of words in the last email
                            num_words, # number of words in the original email
                            gsub("\n|\""," ",temp.event[9])) # content
      }
      interest_record <- rbind(interest_record, evt_record)
    }
  }else{
    temp.event[,4:6] <- tolower(temp.event[,4:6])
    sender_no <- NULL
    receiver_no <- NULL
    for(usr in 1:n.usr){
      temp_id <- identifier_list[[usr]]
      if(all(!is.na(temp.event[1, 4:6]))){
        if(sum(str_detect(temp.event[1, 4],temp_id))>0){
          sender_no <- c(sender_no, usr)
        }
        if(sum(str_detect(temp.event[1, 5],temp_id)) + sum(str_detect(temp.event[1, 6],temp_id)) >0){
          receiver_no <- c(receiver_no, usr)
        }
      }
    }
    if(!is.null(sender_no) & !is.null(receiver_no) & !identical(sender_no,receiver_no)){
      # count the number of words in the email content
      content_last <- temp.event[1,9]
      content_last <- gsub("[^[:alnum:]]", " ", x= content_last)
      num_words_last <- length(unlist(strsplit(x=content_last, split=" +")))
      
      content_ori <- temp.event[n.emails,9]
      content_ori <- gsub("[^[:alnum:]]", " ", x= content_ori)
      num_words_ori <- length(unlist(strsplit(x=content_ori, split=" +")))
      
      if(length(receiver_no) == 1){
        evt_record <- c(evt, # Event index,
                        user[evt], # folder name
                        sender_no, # sender
                        employees[sender_no, 4], # name of sender
                        receiver_no, # receiver
                        employees[receiver_no,4], # name of receiver
                        temp.event[1,1], # Local time
                        temp.event[1,2], # UTC time
                        temp.event[1,3], # Local time
                        length(receiver_no), # Numeber of receiver
                        gsub("\""," ",temp.event[1,7]), # subject
                        temp.event[2,8], # category
                        num_words_last, # number of words in the last email
                        num_words_ori, # number of words in the original email
                        gsub("\n|\""," ",temp.event[1,9])) # content
      }else{
        evt_record <- cbind(evt, # Event index,
                            user[evt], # folder name
                            sender_no, # sender
                            employees[sender_no, 4], # name of sender
                            receiver_no, # receiver
                            employees[receiver_no,4], # name of receiver
                            temp.event[1,1], # Local time
                            temp.event[1,2], # UTC time
                            temp.event[1,3], # Local time
                            length(receiver_no), # Numeber of receiver
                            gsub("\""," ",temp.event[1,7]), # subject
                            temp.event[2,8], # category
                            num_words_last, # number of words in the last email
                            num_words_ori, # number of words in the original email
                            gsub("\n|\""," ",temp.event[1,9])) # content
      }
      interest_record <- rbind(interest_record, evt_record)
    }
  }
  count <- count + 1
  if(evt %% 1000 == 0){
    message(paste0("Finish the first ",count," events."))
    message(paste0("There are ",nrow(interest_record)," interesting events."))
  }
}

# write.csv(interest_record,file = paste0("interested_events/events_among_employees.csv"))
colnames(interest_record) <- c("EventIndex","FolderName","SenderID","SenderName",
                               "ReceiverID","ReceiverName","OriginalTime","UTCtime",
                               "LocalTime","NumReciever","Subject","Category",
                               "LastWordNum","OriginalWordNum","Content")

##############################################################################################
# Obtain the adjacent matrix from interest_record and output the email content for each pair #
##############################################################################################

pairs_in_events <- paste0(interest_record[,3],"-",interest_record[,4])

factor_sender <- factor(interest_record[,3],levels = 1:n.usr)
factor_receiver <- factor(interest_record[,4],levels = 1:n.usr)

n.interestevents <- nrow(interest_record)

if(!dir.exists("interested_events")){
  dir.create("interested_events")
}


# write out the emails between a pair of users into an individual files
for(evt in 1:n.interestevents){#n.interestevents){
  temp.event <- interest_record[evt,]
  sender <- as.numeric(temp.event[3])
  receiver <- as.numeric(temp.event[5])
  if(sender > receiver){
    fil_name <- paste0("interested_events/",receiver,"-",sender,".txt")
    write.table(matrix(temp.event,nrow = 1), file = fil_name,append = TRUE, 
                row.names = FALSE, sep = ",", col.names = FALSE)
  }else if(receiver > sender){
    fil_name <- paste0("interested_events/",sender,"-",receiver,".txt")
    write.table(matrix(temp.event,nrow = 1), file = fil_name,append = TRUE, 
                row.names = FALSE, sep = ",", col.names = FALSE)
  }
}


############################################################
# Working on each pair of nodes to remove redundent emails #
############################################################
dir_interest<-"interested_events/"
if(!dir.exists(dir_interest)){
  dir.create(dir_interest)
}

setwd(dir_interest)
name.file<-list.files(dir_interest) 
n.file<-length(name.file)
NameNoExtension <- sub(".txt","",name.file)
adj_matrix<- matrix("false",n.usr,n.usr)

#########################################################
# Collect the response or censoring time and covariates #
#########################################################
sur_collect <- NULL

# The threshold to regard an email as unreplied
thres_day <- 21

# sr_time <- NULL
# pair_info<- NULL
# data.i <- read.table(file = name.file[1],header = TRUE)
for(fl in 1:n.file){
  pair<-as.numeric(unlist(strsplit(NameNoExtension[fl],"-")))#i,j
  if(pair[1]!=pair[2]){
    
    temp.events<-read.csv(file = name.file[fl],stringsAsFactors = FALSE,header = FALSE)
    
    colnames(temp.events) <- c("EventIndex","FolderName","SenderID","SenderName",
                               "ReceiverID","ReceiverName","OriginalTime","UTCtime",
                               "LocalTime","NumReciever","Subject","Category",
                               "LastWordNum","OriginalWordNum","Content")
    
    # exclude repeat emails with the same sending time
    unique.events <- temp.events[!duplicated(temp.events[,8]),]
    
    # exclude repeat emails with the same content within 12 hours
    unique.events[which(is.na(unique.events[,11])),11] <- ""
    unique.events[which(is.na(unique.events[,13])),13] <- ""
    
    
    # some repeated emails have contents only dfferent in the number of blanks
    content_no_blank <- gsub(" ","",unique.events[,13])
    
    
    duplicated_content <- which(duplicated(content_no_blank))
    duplicated_index <- NULL
    for(evt in duplicated_content){
      # compare the sending time between the current email and the last same-content email
      index_last <- max(which(content_no_blank[1:(evt-1)]==content_no_blank[evt]))
      if(difftime(unique.events[evt,8],unique.events[index_last,8],units = "hours") < 13){
        duplicated_index <- c(duplicated_index, evt)
      }
    }
    if(!is.null(duplicated_index)){
      unique.events <- unique.events[-duplicated_index, ]
    }
    
    # find the reply time
    # Get the prefix of email subjects
    loc_colon <- unlist(regexec(":",unique.events[,11]))
    num_eml <- nrow(unique.events)
    
    # subject without prefix
    sub_no_prefix <- rep(NA, num_eml)
    sub_prefix <- rep(NA,num_eml)
    
    for(eml in 1:num_eml){
      if(loc_colon[eml] > 0){
        sub_no_prefix[eml] <- substring(unique.events[eml,11],first = loc_colon[eml] + 1)
        sub_prefix[eml] <- substr(unique.events[eml,11], start = 1, stop = loc_colon[eml])
        # remove blanks
        sub_no_prefix[eml] <- gsub(" ","",sub_no_prefix[eml])
        #reply_indicator[eml] <- TRUE
      }else{
        sub_no_prefix[eml] <- unique.events[eml,11]
        sub_prefix[eml] <- ""
        # remove blanks
        sub_no_prefix[eml] <- gsub(" ","",sub_no_prefix[eml])
        #reply_indicator[eml] <- FALSE 
      }
    }
    reply_indicator <- regexpr("re",tolower(sub_prefix)) > 0
    fw_indicator <- regexpr("fw",tolower(sub_prefix)) > 0 | unique.events[,12]=="Forwarded"
    
    #observed survival outcome
    # email indices except the last two emails
    emails_except_2last <- 1:num_eml
    last1 <- NULL
    tmp <- which(as.numeric(unique.events[,3]) == pair[1] & as.numeric(unique.events[,5]) == pair[2])
    if(length(tmp) > 0){
      last1 <- max(tmp)
    }
    last2 <- NULL
    tmp <- which(as.numeric(unique.events[,3]) == pair[2] & as.numeric(unique.events[,5]) == pair[1])
    if(length(tmp) > 0){
      last2 <- max(tmp)
    }
    emails_except_2last <- emails_except_2last[-c(last1,last2)]
    
    # mark original emails as replied
    replied_ind <- rep(NA, num_eml)
    for(eml in which(reply_indicator)){
      # all emails within 21 days before this email can be the original email 
      sender <- as.numeric(unique.events[eml,3])
      receiver <- as.numeric(unique.events[eml,5])
      
      if(eml > 1){
        possible_replied <- which(unique.events[1:(eml - 1),3] == receiver)
        
        if(length(which(sub_no_prefix[possible_replied] == sub_no_prefix[eml])) > 0){
          # consider the closest email with the same subject
          eml_replied <- possible_replied[max(which(sub_no_prefix[possible_replied] == sub_no_prefix[eml]))]
          if(difftime(unique.events[eml,8], unique.events[eml_replied,8], units = "hours") < thres_day * 24){
            replied_ind[eml_replied] <- eml
          }
        }
      }
    }
    
    if(length(emails_except_2last) > 0){
      sur_out <- matrix(0, length(emails_except_2last), 10)
      k1 <- 0
      k2 <- 0
      
      #covariate matrix
      sur_out[,5] <- 1
      sur_out[,6] <- weekdays(as.Date(unique.events[emails_except_2last,7],format = "%a, %d %b %Y %H:%M:%S")) %in% c("Saturday","Sunday")
      sur_out[,7] <- fw_indicator[emails_except_2last]
      sur_out[,8] <- log(as.numeric(unique.events[emails_except_2last,10]))
      sur_out[,9] <- log1p(as.numeric(unique.events[emails_except_2last,14]))
      
      colnames(sur_out) <- c("i","j","k","nu","X0","Weekend","Foward","LogNumReceiver","LogNumWordsPlus1OriginalEmail","Y")
      
      index <- 1
      for(eml in emails_except_2last){
        sender <- as.numeric(unique.events[eml,3])
        receiver <- as.numeric(unique.events[eml,5])
        sur_out[index,1:2] <- c(sender,receiver)
        if(sender > receiver){
          k1 <- k1 + 1
          sur_out[index,3] <- k1
        }else{
          k2 <- k2 + 1
          sur_out[index,3] <- k2
        }
        
        # if it is a replied email, return a reply time. Or, return a censoring time as 21 days
        if(!is.na(replied_ind[eml])){
          reply_eml <- replied_ind[eml]
          sur_out[index,4] <- 1
          sur_out[index,10] <- difftime(unique.events[reply_eml,8], unique.events[eml,8], units = "hours")
        }else if(length(which(sub_no_prefix[(eml + 1):num_eml] == sub_no_prefix[eml])) > 0){
          sur_out[index,4] <- 0
          recent_eml <- min(which(sub_no_prefix[(eml + 1):num_eml] == sub_no_prefix[eml]))
          sur_out[index,10] <- min(thres_day * 24, difftime(unique.events[eml + recent_eml,8], unique.events[eml,8], units = "hours"))
        }else{
          sur_out[index,4] <- 0
          sur_out[index,10] <- thres_day * 24
        }
        
        index <- index + 1
      }
      fil_name <- paste0(NameNoExtension[fl],"_replytime.txt")
      sur_collect <- rbind(sur_collect, sur_out)
      write.csv(sur_out, file = fil_name, row.names = FALSE)
    }
    # save the unique events into file
    fil_name <- paste0(NameNoExtension[fl],"_unique.txt")
    write.table(unique.events, file = fil_name, row.names = FALSE, sep = ",")
  }
  message(paste0("Finish the pair ",NameNoExtension[fl]))
}

# sort reply emails by the user order and include pairs without failure time
tol_rep <- nrow(sur_collect)
p.cov <- ncol(sur_collect) - 6
usr_no_list <- as.numeric(names(table(sur_collect[,1:2])))
N <-length(usr_no_list)

# order the response time by senders and receivers

Dat <- sur_collect
ord_rep <- order(Dat[,1],Dat[,2])
Dat <- Dat[ord_rep,]

# exclude user who do not reply any email and are not replied by others
num.pair <- cbind(num.pair, 0)
index_rep <- 0
for(i in 1:tol_pair){
  nij <- num.pair[i,3]
  num.pair[i,4] <- sum(Dat[index_rep + 1:nij, 4])
  index_rep <- index_rep + nij
}

num_fail_usr <- rep(0,max(usr_no_list))
for(s in usr_no_list){
  send_pair <- which(num.pair[,1] == s)
  num_fail_usr[s] <- num_fail_usr[s] + sum(num.pair[send_pair,4])

  rece_pair <- which(num.pair[,2] == s)
  num_fail_usr[s] <- num_fail_usr[s] + sum(num.pair[rece_pair,4])
}

usr_fail <- which(num_fail_usr > 0)

setwd("../")
save(sur_collect,file="Reponse_time_collection.RData")

######################################################################
# Generate the adjacent matrix of the undirected communication graph #
######################################################################
com_pairs <- matrix(as.numeric(unlist(strsplit(NameNoExtension,"-"))),nrow = 2)
n_pairs <- ncol(com_pairs)
com_pairs <- matrix(as.numeric(factor(com_pairs)), nrow = 2)
for(i in 1:n_pairs){
  sender <- com_pairs[1,i]
  receiver <- com_pairs[2,i]
  adj_matrix[sender,receiver] <- "true"
}

dat_numpair_adj <- data.frame(Sender = com_pairs[1,],
                              Receiver = com_pairs[2,])

write.table(adj_matrix, file = "adj_enron_communications.txt", row.names = FALSE, col.names = FALSE, sep = ",", quote = FALSE)


