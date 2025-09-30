setwd("/Users/sxy/Edinburgh/Extended Statstical Programming/Assignment 1")
a <- scan("shakespeare.txt",what="character",skip=83,nlines=196043-83,
          fileEncoding="UTF-8")

##------------Pre-processing-----------------
#-------------Remove Stage Directions--------------
# Find words containing '['
open_brackets <- grep("\\[", a)
a[open_brackets]
remove <- numeric()
for(i in open_brackets){
  #find the words with ']' in next 100 words
  j = i + grep("\\]",a[i:(i+100)])[1] - 1 #since that need to modify the correct position of ']'
  if(!is.na(j)){
    remove = c(remove,i:j)
  }
}


#remove the stage direction
a <- a[-unique(remove)]



#-------------Remove Uppercase Words and Numerals--------------
filter_word <- function(word_vector){
  uppercase <- (toupper(word_vector) == word_vector) &  #Check if a word consists entirely of uppercase letters
       (word_vector != 'I') & (word_vector != 'A') #Except 'I' and 'A'
  numerals <- grepl("^[0-9]+$",word_vector) # Use regular expression to remove the arabic numerals
  return(word_vector[!(uppercase | numerals)])
}

a <- filter_word(a)


#-------------Remove Apostrophes and Hyphens--------------
?gsub

a_original <- a
a <- gsub("_", "", a)
a <- gsub("-", "", a)

identical(a, a_original) #To make sure there exists changes after doing "gsub"
setdiff(a,a_original)

#-------------Create split_punct Function--------------
split_punct <- function(words, punct){
  punct_pattern <- "[,.!?;:]"
  # Directly find words containing punctuation
  punct_words <- grepl(punct_pattern, words)
  ii <- which(punct_words) #Find the position of words with punctuation
  xs <- rep("",length(words)+length(ii))
  i <- 1
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

punct_marks <- c(",", ".", ";", "!", ":", "?")
a <- split_punct(a, punct_marks)

a <- tolower(a)
a

##------------Building vocabulary and word frequency statistics-----------------

unique_words <- unique(a) #find the vector of unique words

index_vector <- match(a, unique_words) 
#Find the index in a corresponding to each word in the unique word vector. The  

word_occurs_time <- tabulate(index_vector)
#Count the number of occurrences of each unique word

ranks <- rank(-word_occurs_time)
#Words with high frequency need to be sorted at the front
b <- unique_words[ranks <= 1000] 
#Extract the top 1000 high-frequency words as vector b



#------------Build the M matrix, token---------------
a_token <- match(a,b) # generate the vector contains the tokens for representing the whole text
create_sequence_matrix <- function(a, tokens, mlag) {
  # check parameters
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
  M <- M[!apply(M, 1, function(row) any(is.na(row))), ]
  return(M)
}

M_4 <- create_sequence_matrix(a, a_token,4)
M <- create_sequence_matrix(a, a_token,4)
#sequence matrix with mlag = 4

#-----------Write the function of next words--------------
next.word = function(key, M, M1, w = rep(1, ncol(M) - 1)){
  # key: Current word sequence (token vector)Current word sequence (token vector)
  # M: Sequence matrix
  # M1: Full text token vector (for degenerate solutions)
  # w: Mixed weight vector (no normalization required)
  
  mlag = ncol(M) - 1  #Model order
  
  # Initialize result container
  u = numeric(0)
  weight = numeric(0)
  for(i in 1:min(length(key),mlag)){
    mc = mlag - i + 1
    mc_to_mlag = mc : mlag
    current_key = key[(length(key) - i + 1):length(key)]
    # Find matching rows - strict NA matching version
    if(any(is.na(current_key))) {
      # Case involving NA: strict matching
      
      if(all(is.na(current_key))) {
        # Scenario 1: All positions are NA → Match rows where all columns are NA
        na_mask <- apply(M[, mc_to_mlag, drop = FALSE], 1, function(x) all(is.na(x)))
        matching_rows <- which(na_mask)
        
      } else {
        # Scenario 1: All positions are NA → Match rows where all columns are NA
        matching <- rep(TRUE, nrow(M))  # Initialize as all matches
        
        for(j in 1:length(current_key)) {
          col_data <- M[, mc_to_mlag[j]]
          
          if(is.na(current_key[j])) {
            # Where the key is NA： The matching token vector must also be NA
            matching <- matching & is.na(col_data)
          } else {
            # All entries in the key are non NA: must match exactly and cannot be NA
            matching <- matching & (col_data == current_key[j]) & !is.na(col_data)
          }
        }
        matching_rows <- which(matching)
      }
      
    } else if(length(current_key) == 1) {
      # Single key matching (without NA)
      matches = M[, mc] == current_key
      matching_rows = which(matches & is.finite(M[, mc]))
      
    } else {
      # Normal multi key matching (no NA)
      ii = colSums(!t(M[, mc_to_mlag, drop = FALSE]) == current_key)
      matching_rows = which(ii == 0 & is.finite(ii))
    }
    
    if (length(matching_rows) > 0) {
      # Extract the next word (mlag+1 column) and remove NA
      next_word = M[matching_rows, mlag + 1]
      next_word = next_word[!is.na(next_word) & is.finite(next_word)]
      
      if (length(next_word) > 0) {
        # Assign weights to each matching word
        # Weight=w [i] * (1/number of context matches)
        n_matches = length(next_word)
        
        # Collect words and weights
        u <- c(u, next_word)
        weight <- c(weight, rep(w[i] / n_matches, n_matches))
      }
    }
  }
  # Dealing with degraded solutions: If no match is found
  if (length(u) == 0) {
    # Randomly select from common words throughout the text (based on word frequency)
    valid_tokens <- M1[is.finite(M1) & !is.na(M1)]
    if (length(valid_tokens) > 0) {
      token_counts <- table(valid_tokens)
      common_tokens <- as.numeric(names(token_counts))
      freq_probs <- as.numeric(token_counts) / sum(token_counts)
      
      return(sample(common_tokens, 1, prob = freq_probs))
    } else {
      # Extreme case: Randomly return a valid token
      valid_indices <- which(is.finite(M1) & !is.na(M1))
      if (length(valid_indices) > 0) {
        return(sample(M1[valid_indices], 1))
      } else {
        return(1)  # Default return of the first word
      }
    }
  }
  
  # Merge weights of the same token in u (add directly)
  unique_u = unique(u)
  combined_weights = sapply(unique_u, function(token) {
    sum(weight[u == token])
  })
  
  # Note: Normalized probability is required for sampling here, but the weight w itself does not need to be normalized
  # This is a necessary step before sampling, not the normalization of w
  if (sum(combined_weights) > 0) {
    probs = combined_weights / sum(combined_weights)
  } else {
    probs = rep(1/length(unique_u), length(unique_u))
  }
  
  return(sample(unique_u, 1, prob = probs))
}



# ----------------------Q8:Choose a single starting word (non punctuation)----------------
select_start_token <- function(b, M1, specific_word = NULL) {
  # b: Vocabulary vector
  # M1: Full text token vector
  # specific_word: Optional, specify the starting word such as' romeo '
  
  if (!is.null(specific_word)) {
    # Use the specified starting word
    specific_word <- tolower(specific_word)
    start_token <- which(b == specific_word)
    
    if (length(start_token) > 0) {
      cat("Use specified starting words: '", specific_word, "' (token ", start_token, ")\n", sep = "")
      return(start_token)
    } else {
      cat("'", specific_word, "'not in the vocabulary, use random selection\n", sep = "")
    }
  }
  
  # Randomly select starting words (excluding punctuation marks)
  valid_tokens <- M1[!is.na(M1) & is.finite(M1)]
  
  # Exclude punctuation marks
  punctuation <- c(".", ",", ";", "!", "?", ":", "-", "--", "(", ")", "[", "]")
  non_punct_tokens <- valid_tokens[!b[valid_tokens] %in% punctuation]
  
  if (length(non_punct_tokens) > 0) {
    start_token <- sample(non_punct_tokens, 1)
    cat("Randomly select starting words: '", b[start_token], "' (token ", start_token, ")\n", sep = "")
    return(start_token)
  } else {
    # If there are no non punctuation words, use any valid words
    start_token <- sample(valid_tokens, 1)
    cat("Randomly select starting words: '", b[start_token], "' (token ", start_token, ")\n", sep = "")
    return(start_token)
  }
}

# Test for Question 8
cat("=== Test for Question 8 ===\n")
start_token1 <- select_start_token(b, a_token)  # randomly select
start_token2 <- select_start_token(b, a_token, "romeo")  # fixed word
start_token3 <- select_start_token(b, a_token, "love")   # fixed word



#--------------------- Q9: Generate a complete sentence until encountering a period----------------------
generate_sentence_until_fullstop <- function(start_token, M, M1, b, w = rep(1, ncol(M) - 1), max_words = 100) {
  # Start_token: Starting word token
  # Other parameters are the same as the next. word function
  
  mlag <- ncol(M) - 1
  
  # Initialization
  current_sequence <- start_token
  generated_tokens <- start_token
  word_count <- 1
  
  cat("Start generating sentences...\n")
  cat("Starting words: '", b[start_token], "'\n", sep = "")
  cat("Generating process:\n", b[start_token])
  
  # Generate until a period is encountered or the maximum length is reached
  for (i in 1:max_words) {
    # Generate the next word using the next.word function
    next_token <- next.word(current_sequence, M, M1, w)
    generated_tokens <- c(generated_tokens, next_token)
    word_count <- word_count + 1
    
    # Display the generating process
    cat(" ", b[next_token])
    
    # Update the current sequence (keeping the length within mlag)
    if (length(current_sequence) >= mlag) {
      current_sequence <- c(current_sequence[2:mlag], next_token)
    } else {
      current_sequence <- c(current_sequence, next_token)
    }
    
    # Check if a period has been encountered
    if (b[next_token] == ".") {
      cat(" [Encounter period, end]\n")
      break
    }
    
    # Display line breaks every 10 words
    if (i %% 10 == 0) cat("\n")
  }
  
  # Check if the maximum length has been reached but no period has been encountered
  if (b[generated_tokens[length(generated_tokens)]] != ".") {
    cat(" [Reached maximum length, no period encountered]\n")
  }
  
  # Convert the token vector back to words and beautify the output
  sentence_words <- b[generated_tokens]
  sentence <- paste(sentence_words, collapse = " ")
  
  # Processing punctuation spaces (removing spaces before punctuation)
  sentence <- gsub("\\s+([,.!?;:])", "\\1", sentence)
  
  # Capitalize the first letter
  sentence <- paste0(toupper(substr(sentence, 1, 1)), substr(sentence, 2, nchar(sentence)))
  
  # Ensure to end with a period (if not already)
  if (!grepl("[.!?]$", sentence)) {
    sentence <- paste0(sentence, ".")
  }
  
  cat("Generation completed, ", word_count, "words in total.\n")
  return(list(
    sentence = sentence,
    tokens = generated_tokens,
    word_count = word_count
  ))
}

# Complete testing of Question 8 and Question 9
cat("=== Complete testing for questions 8 and 9 ===\n\n")


cat("Question 8: Choose the starting word\n")


cat("Method 1: Random selection:\n")
start_random <- select_start_token(b, a_token)

#Method 2: Specify interesting vocabulary
cat("\nMethod 2- Specify Vocabulary:\n")
start_romeo <- select_start_token(b, a_token, "romeo")
start_love <- select_start_token(b, a_token, "love")
start_king <- select_start_token(b, a_token, "king")

# Question 9: Generate sentences
cat("\nQuestion 9: Generate complete sentences\n")

cat("1. Markov Model Generation:\n")


sentences_markov <- list()

cat("Sentence 1 (random start):\n")
sentence1 <- generate_sentence_until_fullstop(start_random, M, a_token, b)
cat("Generated sentence:\n", sentence1$sentence, "\n\n")
sentences_markov[["Random start"]] <- sentence1$sentence

cat("Sentence2 (starting word 'romeo'):\n")
sentence2 <- generate_sentence_until_fullstop(start_romeo, M, a_token, b)
cat("Generated sentence", sentence2$sentence, "\n\n")
sentences_markov[["starting word 'romeo'"]] <- sentence2$sentence

cat("Sentence 3 (starting word 'love'):\n")
sentence3 <- generate_sentence_until_fullstop(start_love, M, a_token, b)
cat("Generated sentence:\n", sentence3$sentence, "\n\n")
sentences_markov[["Starting word 'love'"]] <- sentence3$sentence


