# read icd10 and birth file ----------------------------------------------------

df <- data.table::fread("./data/icd10withtime_extract.csv")

# read link file ---------------------------------------------------------------

link <- data.table::fread("./data/linker_app15825_withexcl.csv")
colnames(link) <- c("ieu","eid")

df_icd10 = df[,c('eid',paste0("41270-0.",seq(0,225)))]

codes <- reshape2::melt(df_icd10,
                        id.vars = c("eid"),
                        variable.name = "var_code",
                        value.name = "value_code")
codes$var_code <- gsub("-.*","",as.character(codes$var_code))

# find id for hypertension (I10,I15) for case ----------------------------------

codes_i1015 = filter(codes,value_code%in%c('I10',"I150","I151","I152","I158","I159"))
idhpt_cases = sort(unique(codes_i1015$eid))

df_hpt = data.frame(df[,"eid"])
df_hpt$hpt = ifelse(df_hpt$eid %in%idhpt_cases,1,0)
df_hpt = merge(df_hpt,link,by =c("eid"))

df_hpt = df_hpt[,c(3,3,2)]
colnames(df_hpt) = c("FID","IID","hpt")

# write hypertension input file ------------------------------------------------

write.table(df_hpt,"./data/df_hpt.txt",
            sep = " ",
            row.names = F, col.names = T,quote =F)

# for early- and late- onset hypertension input --------------------------------
# read data --------------------------------------------------------------------

df <- data.table::fread("./data/icd10withtime_extract.csv")

# read link file ---------------------------------------------------------------

link <- data.table::fread("./data/linker_app15825_withexcl.csv")
colnames(link) <- c("ieu","eid")

df_icd10 = df[,c('eid',paste0("41270-0.",seq(0,225)))]

df_icd10date = df[,c('eid',paste0("41280-0.",seq(0,225)))]

df_birth = df[,c("eid","34-0.0", "52-0.0")]

# reshape the data -------------------------------------------------------------

codes <- reshape2::melt(df_icd10,
                        id.vars = c("eid"),
                        variable.name = "var_code",
                        value.name = "value_code")

codes$var_code <- gsub("-.*","",as.character(codes$var_code))

dates <- reshape2::melt(df_icd10date,
                        id.vars = c("eid"),
                        variable.name = "var_date",
                        value.name = "value_date")

dates$var_date <- gsub("-.*","",as.character(dates$var_date))

code_date = cbind(codes,dates)

code_date_case = subset(code_date, value_code%in%c('I10',"I150","I151","I152","I158","I159"))

# randomly give the day of birth for each case ---------------------------------
# identify the leap year since 1900 --------------------------------------------

colnames(df_birth) =c("eid","Year","Month")
leap_year = c('2020','2016','2012','2008',
              '2004','2000','1996','1992','1988',
              '1984','1980','1976','1972','1968',
              '1964','1960','1956','1952','1948',
              '1944','1940','1936','1932','1928',
              '1924','1920','1916','1912','1908',
              '1904')
df_ly = subset(df_birth, Year %in%leap_year)
df_nly = subset(df_birth, Year %notin%leap_year)

month_g1 = c('1','3','5','7','8','10','12')
month_g2 = c('4','6','9','11')

set.seed(321)
df_ly_g1 = subset(df_ly, Month%in%month_g1)
df_ly_g1$Day = sample(1:31,75493,replace = T)
df_ly_g2 = subset(df_ly, Month%in%month_g2)
df_ly_g2$Day = sample(1:30,41793,replace = T)
df_ly_2 = subset(df_ly, Month == 2)
df_ly_2$Day = sample(1:29,10351,replace = T)
df_nly_g1 = subset(df_nly, Month%in%month_g1)
df_nly_g1$Day = sample(1:31,222213,replace = T)
df_nly_g2 = subset(df_nly, Month%in%month_g2)
df_nly_g2$Day = sample(1:30,122981,replace = T)
df_nly_2 = subset(df_nly, Month == 2)
df_nly_2$Day = sample(1:28,29590,replace = T)

df_birth = rbind(df_ly_g1,df_ly_g2,df_ly_2,
                 df_nly_g1,df_nly_g2,df_nly_2)

# assign birth date to cases ---------------------------------------------------

df_birth$date = make_date(year = df_birth$Year, month = df_birth$Month, day = df_birth$Day)
df_birth_I10date = left_join(code_date_case,df_birth,by = 'eid')
df_birth_I10date$value_date = as.Date(df_birth_I10date$value_date)

# calculate the difference between birth and 1st diagnosis ---------------------

df_birth_I10date$diff = as.numeric(df_birth_I10date$value_date - df_birth_I10date$date)/365.24219

# age threshold for 50 - 60
df_birth_I10date$early50 <- ifelse(df_birth_I10date$diff <=50,1,0)
df_birth_I10date$late60 <- ifelse(df_birth_I10date$diff >60,1,0)

df_birth_I10date = merge(df_birth_I10date,link,by =c("eid"))

early50id = subset(df_birth_I10date,early50 ==1)$ieu

late60id =  subset(df_birth_I10date,late60 ==1)$ieu

# read the whole sample file ---------------------------------------------------

df_hpt = data.table::fread("./data/df_hpt.txt")
df_hpt$early50 = ifelse(df_hpt$IID%in%early50id,1,0)
df_hpt$late60 = ifelse(df_hpt$IID%in%late60id,1,0)

df_early50 = subset(df_hpt, hpt == 0 |early50 == 1)[,c("FID","IID","early50")]

df_late60 = subset(df_hpt, hpt == 0 |late60 == 1)[,c("FID","IID","late60")]

# write input files for early- and late- onset hypertension  -------------------

write.table(df_early50,"./data/df_early50.txt",
            sep = " ",
            row.names = F, col.names = T,quote =F)

write.table(df_late60,"./data/df_late60.txt",
            sep = " ",
            row.names = F, col.names = T,quote =F)
