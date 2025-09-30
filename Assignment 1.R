#setwd("/Users/sxy/Edinburgh/Extended Statstical Programming/Assignment 1")
a <- scan("shakespeare.txt",what="character",skip=83,nlines=196043-83,
          fileEncoding="UTF-8")

##------------Pre-processing-----------------
#-------------Remove Stage Directions--------------
# Find words containing '['
open_brackets = grep("\\[", a)
a[open_brackets]
remove = unmatched_brackets = numeric()
for(i in open_brackets){
  #find the words with ']' in next 100 words
  j = i + grep("\\]",a[i:(i+100)])[1] - 1 #since that need to modify the correct position of ']'
  if(!is.na(j)){
    remove = c(remove,i:j)
  }else {                  #if j is NA, which means there is no "]" after "[" in 100 characters
    unmatched_brackets = c(unmatched_brackets, i)
    cat("The position of unmatched brackets:", i, "The text is:", a[i], "\n") 
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
a = gsub("_", "", a)
a = gsub("-", "", a)

identical(a, a_original) #To make sure there exists changes after doing "gsub"
setdiff(a,a_original)

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
a_token = match(a,b) # generate the vector contains the tokens for representing the whole text
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
  M <- M[!apply(M, 1, function(row) any(is.na(row))), ] #remove NA in the first and last row
  return(M)
}

M_4 = create_sequence_matrix(a, a_token,4)
M_4 #sequence matirx with mlag = 4

#-----------Write the function of next words--------------
next.word = function(key, M, M1, w = rep(1, ncol(M) - 1)){
  # key: 当前词序列（token向量）
  # M: 序列矩阵
  # M1: 全文token向量（用于退化解）
  # w: 混合权重向量（不需要归一化）
  
  mlag = ncol(M) - 1  # 模型阶数
  
  # 初始化结果容器
  u = numeric(0)
  weight = numeric(0)
  for(i in 1:min(length(key),mlag)){
    mc = mlag - i + 1
    mc_to_mlag = mc : mlag
    current_key = key[(length(key) - i + 1):length(key)]
    ## 找匹配的行 - 严格NA匹配版本
    if(any(is.na(current_key))) {
      # 包含NA的情况：严格匹配
      
      if(all(is.na(current_key))) {
        # 情况1：所有位置都是NA → 匹配所有列都是NA的行
        na_mask <- apply(M[, mc_to_mlag, drop = FALSE], 1, function(x) all(is.na(x)))
        matching_rows <- which(na_mask)
        
      } else {
        # 情况2：部分位置是NA → 严格匹配
        matching <- rep(TRUE, nrow(M))  # 初始化为全部匹配
        
        for(j in 1:length(current_key)) {
          col_data <- M[, mc_to_mlag[j]]
          
          if(is.na(current_key[j])) {
            # key中是NA：数据中也必须是NA
            matching <- matching & is.na(col_data)
          } else {
            # key中是非NA：必须精确匹配，且不能是NA
            matching <- matching & (col_data == current_key[j]) & !is.na(col_data)
          }
        }
        matching_rows <- which(matching)
      }
      
    } else if(length(current_key) == 1) {
      # 单键匹配（无NA）
      matches = M[, mc] == current_key
      matching_rows = which(matches & is.finite(M[, mc]))
      
    } else {
      # 正常多键匹配（无NA）
      ii = colSums(!t(M[, mc_to_mlag, drop = FALSE]) == current_key)
      matching_rows = which(ii == 0 & is.finite(ii))
    }
    
    if (length(matching_rows) > 0) {
      # 提取下一个词（第mlag+1列），去除NA
      next_word = M[matching_rows, mlag + 1]
      next_word = next_word[!is.na(next_word) & is.finite(next_word)]
      
      if (length(next_word) > 0) {
        # 为每个匹配的word分配权重
        # 权重 = w[i] * (1/该上下文匹配次数)
        n_matches = length(next_word)
        
        # 收集words和权重
        u <- c(u, next_word)
        weight <- c(weight, rep(w[i] / n_matches, n_matches))
      }
    }
  }
  # 处理退化解：如果没有找到任何匹配
  if (length(u) == 0) {
    # 从全文常见词中随机选择（基于词频）
    valid_tokens <- M1[is.finite(M1) & !is.na(M1)]
    if (length(valid_tokens) > 0) {
      token_counts <- table(valid_tokens)
      common_tokens <- as.numeric(names(token_counts))
      freq_probs <- as.numeric(token_counts) / sum(token_counts)
      
      return(sample(common_tokens, 1, prob = freq_probs))
    } else {
      # 极端情况：随机返回一个有效token
      valid_indices <- which(is.finite(M1) & !is.na(M1))
      if (length(valid_indices) > 0) {
        return(sample(M1[valid_indices], 1))
      } else {
        return(1)  # 默认返回第一个词
      }
    }
  }
  
  # 合并相同token的权重（直接相加）
  unique_u = unique(u)
  combined_weights = sapply(unique_u, function(token) {
    sum(weight[u == token])
  })
  
  # 注意：这里需要归一化概率用于抽样，但权重w本身不需要归一化
  # 这是抽样前的必要步骤，不是对w的归一化
  if (sum(combined_weights) > 0) {
    probs = combined_weights / sum(combined_weights)
  } else {
    probs = rep(1/length(unique_u), length(unique_u))
  }
  
  return(sample(unique_u, 1, prob = probs))
}



# ----------------------第8问：选择单个起始词（非标点）----------------
select_start_token <- function(b, M1, specific_word = NULL) {
  # b: 词汇表向量
  # M1: 全文token向量
  # specific_word: 可选，指定起始词如"romeo"
  
  if (!is.null(specific_word)) {
    # 使用指定的起始词
    specific_word <- tolower(specific_word)
    start_token <- which(b == specific_word)
    
    if (length(start_token) > 0) {
      cat("使用指定起始词: '", specific_word, "' (token ", start_token, ")\n", sep = "")
      return(start_token)
    } else {
      cat("指定词 '", specific_word, "' 不在词汇表中，使用随机选择\n", sep = "")
    }
  }
  
  # 随机选择起始词（排除标点符号）
  valid_tokens <- M1[!is.na(M1) & is.finite(M1)]
  
  # 排除标点符号
  punctuation <- c(".", ",", ";", "!", "?", ":", "-", "--", "(", ")", "[", "]")
  non_punct_tokens <- valid_tokens[!b[valid_tokens] %in% punctuation]
  
  if (length(non_punct_tokens) > 0) {
    start_token <- sample(non_punct_tokens, 1)
    cat("随机选择起始词: '", b[start_token], "' (token ", start_token, ")\n", sep = "")
    return(start_token)
  } else {
    # 如果没有非标点词，使用任何有效词
    start_token <- sample(valid_tokens, 1)
    cat("随机选择起始词: '", b[start_token], "' (token ", start_token, ")\n", sep = "")
    return(start_token)
  }
}

# 测试第8问
cat("=== 第8问测试 ===\n")
start_token1 <- select_start_token(b, a_token)  # 随机选择
start_token2 <- select_start_token(b, a_token, "romeo")  # 指定词
start_token3 <- select_start_token(b, a_token, "love")   # 指定词



#--------------------- 第9问：生成完整句子直到遇到句----------------------
generate_sentence_until_fullstop <- function(start_token, M, M1, b, w = rep(1, ncol(M) - 1), max_words = 100) {
  # start_token: 起始词token
  # 其他参数同next.word函数
  
  mlag <- ncol(M) - 1
  
  # 初始化
  current_sequence <- start_token
  generated_tokens <- start_token
  word_count <- 1
  
  cat("开始生成句子...\n")
  cat("起始词: '", b[start_token], "'\n", sep = "")
  cat("生成过程: ", b[start_token])
  
  # 生成直到遇到句号或达到最大长度
  for (i in 1:max_words) {
    # 使用next.word函数生成下一个词
    next_token <- next.word(current_sequence, M, M1, w)
    generated_tokens <- c(generated_tokens, next_token)
    word_count <- word_count + 1
    
    # 显示生成过程
    cat(" ", b[next_token])
    
    # 更新当前序列（保持长度不超过mlag）
    if (length(current_sequence) >= mlag) {
      current_sequence <- c(current_sequence[2:mlag], next_token)
    } else {
      current_sequence <- c(current_sequence, next_token)
    }
    
    # 检查是否遇到句号
    if (b[next_token] == ".") {
      cat(" [遇到句号，结束]\n")
      break
    }
    
    # 每10个词换行显示
    if (i %% 10 == 0) cat("\n")
  }
  
  # 检查是否达到最大长度但未遇到句号
  if (b[generated_tokens[length(generated_tokens)]] != ".") {
    cat(" [达到最大长度，未遇到句号]\n")
  }
  
  # 将token转换回单词并美化输出
  sentence_words <- b[generated_tokens]
  sentence <- paste(sentence_words, collapse = " ")
  
  # 处理标点符号的空格（移除标点前的空格）
  sentence <- gsub("\\s+([,.!?;:])", "\\1", sentence)
  
  # 首字母大写
  sentence <- paste0(toupper(substr(sentence, 1, 1)), substr(sentence, 2, nchar(sentence)))
  
  # 确保以句号结束（如果还没有）
  if (!grepl("[.!?]$", sentence)) {
    sentence <- paste0(sentence, ".")
  }
  
  cat("生成完成，总共", word_count, "个词。\n")
  return(list(
    sentence = sentence,
    tokens = generated_tokens,
    word_count = word_count
  ))
}

# 完整的第8问和第9问测试
cat("=== 第8问和第9问完整测试 ===\n\n")

# 第8问：选择起始词
cat("第8问：选择起始词\n")
cat(rep("-", 50), "\n")

# 方法1：随机选择
cat("方法1 - 随机选择:\n")
start_random <- select_start_token(b, a_token)

# 方法2：指定有趣词汇
cat("\n方法2 - 指定词汇:\n")
start_romeo <- select_start_token(b, a_token, "romeo")
start_love <- select_start_token(b, a_token, "love")
start_king <- select_start_token(b, a_token, "king")

# 第9问：生成句子
cat("\n第9问：生成完整句子\n")
cat(rep("=", 60), "\n")

# 使用马尔可夫模型生成
cat("1. 马尔可夫模型生成:\n")
cat(rep("-", 40), "\n")

sentences_markov <- list()

cat("句子1 (随机起始):\n")
sentence1 <- generate_sentence_until_fullstop(start_random, M, a_token, b)
cat("结果:", sentence1$sentence, "\n\n")
sentences_markov[["随机起始"]] <- sentence1$sentence

cat("句子2 (起始词'romeo'):\n")
sentence2 <- generate_sentence_until_fullstop(start_romeo, M, a_token, b)
cat("结果:", sentence2$sentence, "\n\n")
sentences_markov[["romeo起始"]] <- sentence2$sentence

cat("句子3 (起始词'love'):\n")
sentence3 <- generate_sentence_until_fullstop(start_love, M, a_token, b)
cat("结果:", sentence3$sentence, "\n\n")
sentences_markov[["love起始"]] <- sentence3$sentence

