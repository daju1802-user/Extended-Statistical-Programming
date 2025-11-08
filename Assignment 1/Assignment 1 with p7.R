setwd("/Users/sxy/Edinburgh/Extended Statstical Programming/Assignment 1")
a <- scan("shakespeare.txt",what="character",skip=83,nlines=196043-83,
          fileEncoding="UTF-8")

##------------Pre-processing-----------------
#-------------Remove Stage Directions--------------
# Find words containing '['
open_brackets <- grep("\\[", a)
a[open_brackets]
remove = numeric()
bracket_info =  data.frame(
  start_index = numeric(),
  end_index = numeric(),
  word_count = numeric(),
  stringsAsFactors = FALSE
)
for(i in open_brackets){
  #find the words with ']' in next 100 words
  j = i + grep("\\]",a[i:(i+100)])[1] - 1
  #since that need to modify the correct position of first ']' in the next 100 words
  if(!is.na(j)){
    remove = c(remove,i:j)#continuously record positions of '[', ']' in the loop
    
    word_count <- j - i + 1
    bracket_info <- rbind(bracket_info, data.frame(
      start_index = i,
      end_index = j,
      word_count = word_count
    ))
    #record the number of words in the deleted brackets
  }
}   
#remove the stage direction
bracket_info#show recorded information of brackets
a = a[-unique(remove)]#delete repeated index to avoid error

freq_table <- table(bracket_info$word_count)
mode_count <- as.numeric(names(freq_table)[which.max(freq_table)])
mode_count

# deal with the unmatched brackets
remove_unmatched <- numeric()

# 
remaining_open <- grep("\\[", a)
if(length(remaining_open) > 0){
  open_info <- data.frame(
    bracket_index = remaining_open,
    bracket_type = "[",
    stringsAsFactors = FALSE
  )
  
  for(i in remaining_open){
    # 
    end_pos <- min(i + mode_count - 1, length(a))
    remove_unmatched <- c(remove_unmatched, i:end_pos)
  }
}

# 
remaining_close <- grep("\\]", a)
if(length(remaining_close) > 0){
  close_info <- data.frame(
    bracket_index = remaining_close,
    bracket_type = "]",
    stringsAsFactors = FALSE
  )
  
  for(i in remaining_close){
    # 
    start_pos <- max(i - mode_count + 1, 1)
    remove_unmatched <- c(remove_unmatched, start_pos:i)
  }
}
remove_unmatched
a <- a[-unique(remove_unmatched)]
a


#-------------Remove Uppercase Words and Numerals-------------
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


#--------5--------
unique_words <- unique(a) #find the vector of unique words

index_vector <- match(a, unique_words) 
#Find the index in a corresponding to each word in the unique word vector. The  

word_occurs_time <- tabulate(index_vector)
#Count the number of occurrences of each unique word

ranks <- rank(-word_occurs_time)
#Words with high frequency need to be sorted at the front
b <- unique_words[ranks <= 1000] 
#Extract the top 1000 high-frequency words as vector b


#-----------6-----------
a_token = match(a,b) # generate the vector contains the tokens for representing the whole text
create_sequence_matrix <- function(a, tokens, mlag) {
  # 参数检查
  if (mlag < 1) {
    stop("mlag should be atleast 1")
  }
  
  n <- length(a)
  if (n <= mlag) {
    stop("the text length must be greater than the mlag value")
  }
  
  
  # Create sequence matrix M
  # Initialize matrix M, with a size of (n-mlag) × (mlag+1)
  M <- matrix(NA, nrow = n - mlag, ncol = mlag + 1)
  
  # Fill each column of the matrix
  for (i in 0:mlag) {
    # The i-th column is the token vector that moves i positions to the right
    M[, i + 1] <- tokens[(1 + i):(n - mlag + i)]
  }
  
 
  
  return(M)
}
M_4 = create_sequence_matrix(a, a_token,4)
M_4 #sequence matrix with mlag = 4


#--------7--------
next.word <- function(key, M, M1, w = rep(1, ncol(M) - 1)) {
  # Parameter description:
  # key: Keyword sequence (word tag vector)
  # M: Sequence matrix
  # M1: The token vector of the entire text
  # w: Mixed weight vector, default equal weight
  
  mlag <- ncol(M) - 1  # Maximum lag value
  
  # Processing key length
  if (length(key) > mlag) {
    # If the key is too long, only use the last mlag words
    key <- key[(length(key) - mlag + 1):length(key)]
    mc <- 1  # Match from the first column onwards
  } else {
    # If the key is short, use a reduced order model
    mc <- mlag - length(key) + 1  # Calculate the starting column
  }
  
  # Get the submatrix of M (columns to match)
  M_sub <- M[, mc:mlag, drop = FALSE]
  
  # Calculate matching degree: For each row, calculate the number of mismatches with the key
  # If all columns in a row match the key, ii[j] = 0
  ii <- colSums(!(t(M_sub) == key))
  
  # Find a perfectly matched row (ii=0 and is a finite value)
  matched_rows <- which(ii == 0 & is.finite(ii))
  
  # If no match is found, use a reduced order model
  if (length(matched_rows) == 0) {
    if (length(key) > 1) {
      # Recursive call, using shorter keys (removing the first word)
      return(next.word(key[-1], M, M1, w))
    } else {
      # If the key length is 1 and there is no match, return the most frequent word
      word_counts <- tabulate(M1)
      return(which.max(word_counts))
    }
  }
  
  # Extract the next column of the matching row (predicted word)
  next_words <- M[matched_rows, mlag + 1]
  
  # Remove NA
  next_words <- next_words[!is.na(next_words)]
  
  if (length(next_words) == 0) {
    # If there is no valid next word, return the most frequent word
    word_counts <- tabulate(M1)
    return(which.max(word_counts))
  }
  
  # Use table() instead of tabulate() to obtain the frequency of each unique word
  word_freq <- table(next_words)
  unique_words <- as.numeric(names(word_freq))
  frequencies <- as.vector(word_freq)
  
  # Ensure that the length of the probability vector is consistent with the length of the sampling object
  if (length(unique_words) > 1) {
    # Calculate probability
    probs <- frequencies / sum(frequencies)
    if (length(probs) == length(unique_words)) {
      next_token <- sample(unique_words, 1, prob = probs)
    } else {
      # If the lengths do not match, use uniform probability
      next_token <- sample(unique_words, 1)
    }
  } else {
    # If there is only one candidate word, return directly
    next_token <- unique_words[1]
  }
  
  return(next_token)
}
key=c(111,99,475,177)
next_word = next.word(key,M_4,a_token)
next_word
b[c(key,next_word)]
