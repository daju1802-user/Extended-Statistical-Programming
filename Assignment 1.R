setwd("/Users/sxy/Edinburgh/Extended Statstical Programming/Assignment 1")
a <- scan("shakespeare.txt",what="character",skip=83,nlines=196043-83,
          fileEncoding="UTF-8")

##------------Pre-processing-----------------
#-------------Remove Stage Directions--------------
# Find words containing '['
open_brackets <- grep("\\[", a)
a[open_brackets]
remove = numeric()
for(i in open_brackets){
  #find the words with ']' in next 100 words
  j = i + grep("\\]",a[i:(i+100)])[1] - 1 #since that need to modify the correct position of ']'
  if(!is.na(j)){
    remove = c(remove,i:j)
  }
}
#remove the stage direction
a = a[-unique(remove)]



#-------------Remove Uppercase Words and Numerals--------------
filter_word = function(word_vector){
  uppercase = (toupper(word_vector) == word_vector) &  #Check if a word consists entirely of uppercase letters
       (word_vector != 'I') & (word_vector != 'A') #Except 'I' and 'A'
  numerals = grepl("^[0-9]+$",word_vector) # Use regular expression to remove the arabic numerals
  return(word_vector[!(uppercase | numerals)])
}

a = filter_word(a)


#-------------Remove Apostrophes and Hyphens--------------
?gsub

a_original <- a
a <- gsub("_", "", a)
a <- gsub("-", "", a)

identical(a, a_original) #To make sure there exists changes after doing "gsub"


#-------------Create split_punct Function--------------
split_punct = function(words, punct){
  punct_pattern = "[,.!?;:]"
  # Directly find words containing punctuation
  punct_words = grepl(punct_pattern, words)
  ii = which(punct_words) #Find the position of words with punctuation
  xs = rep("",length(words)+length(ii))
  i = 1
  for(w in words){
    if(grepl(punct_pattern,w)){
      word_part = gsub(punct_pattern,"",w)
      punct_part = gsub("[^,.!?;:]", "", w) # "^" means not. Not ",.!?;:" in here would be removed
      xs[i] = word_part
      xs[i+1] = punct_part
      i = i + 2
    }else{
      xs[i] = w
      i = i + 1
    }
  }
  return(xs)
}

punct_marks = c(",", ".", ";", "!", ":", "?")
a = split_punct(a, punct_marks)

a <- tolower(a)


#--------第5步--------
unique_words <- unique(a) #find the vector of unique words

index_vector <- match(a, unique_words) 
#Find the index in b corresponding to each word
#The frequency can be counted by numerical operation

word_occurs_time <- tabulate(index_vector)
#Count the number of occurrences of each unique word

ranks <- rank(-word_occurs_time)
#Words with high frequency need to be sorted at the front
b_top1000 <- b[ranks <= 1000]
#Extract the top 1000 high-frequency words from b
