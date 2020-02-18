#' Function to parse from a character vector all character matching keywords from a query using TF-IDF metric.
#' 
#' @param query character or character vector.
#' @param db a character vector to parse query from. If remakedb=FALSE it should be a matrix with each extracted word of db as second column and with corresonding full character of db as first column.
#' @param remakedb a boolean value, default is TRUE. Put to FALSE (you should not do this) if you can provide a matrix as db as describe in db parameter.  
#' @return a character vector with all matches for one query.
#' @details 
#' This function will look for matches between a query and a character database using TF-IDF for text mining.
#' This function is not perfect and require some good practices:
#' - Take look at your database to design your query.
#' - Include both litteral and acronym form (if it exist) of a term in distinct query character.
#' - Avoid using plural forms in query.
#' - Correclty spell your query.
#' - Database will not always be clean, and concatenated words (for instance: forinstance)
#' will not be detected but with a specific corresponding concatened query.
#' @export
#' 
ctQueryTFIDF = function(query, db, threshold=1, remakedb=TRUE) {
    
    # Parse words in db if needed
    if (remakedb) {
        new_db = matrix(ncol = 2)
        for (i in db) {
            ibis = gsub('([[:upper:]][[:lower:]])', ' \\1', i)
            words = as.vector(str_split(ibis, regex("[[:punct:] ]+"), simplify = TRUE))
            one_ct = matrix(data = cbind(rep(i,length(words)), words), ncol=2, nrow = length(words))
            new_db = rbind(new_db, one_ct)
        }
        colnames(new_db) = c("book", "word")
        new_db = as.data.frame(new_db[-c(1,which(new_db[,2] == "")),])
        
        # Homogenize word spelling:
        for (wd in unique(new_db[,2])) {
            word.reg = paste0("^", wd, "s?$")
            db.match.word = str_detect(new_db[,2], regex(word.reg, ignore_case = TRUE))
            # May need optimazation
            new_db[which(db.match.word),2] = wd
        }
        
        # Compute word count and tf_idf
        new_db = new_db %>% count(book, word)
        total_words = new_db %>% 
            group_by(book) %>% 
            summarize(total = sum(n))
        new_db = left_join(new_db, total_words)
        new_db = new_db %>%
            bind_tf_idf(word, book, n)    
        
        # Create books x words tf_idf matrix
        db.words = unique(new_db$word)
        db.names = unique(new_db$book)
        tfidf.table = matrix(0L, nrow = length(db.words), ncol = length(db.names))
        colnames(tfidf.table) = sort(db.names)
        rownames(tfidf.table) = sort(db.words)
        for (i in 1:nrow(new_db)) {
            tfidf.table[as.vector(new_db$word[i]), as.vector(new_db$book[i])] = new_db$tf_idf[i]
        }
        
    } else {
        tfidf.table = db
    }
    # If user provides several queries for one db, apply this function for each query
    if (length(query) > 1) {
        db.match.list = apply(as.array(query), MARGIN = 1, FUN = "ctQueryTFIDF",
                              db=tfidf.table, threshold=threshold, remakedb=FALSE)
        db.match = unlist(db.match.list)
        if (length(which(duplicated(db.match))) > 0) {
            db.match = db.match[-which(duplicated(db.match))]
        }
    } else {
        
        # Parse words in query and find match in db for each words
        query.words = gsub('([[:upper:]][[:lower:]])', ' \\1', query)
        query.words = as.vector(str_split(query.words, regex("[[:punct:] ]+"), simplify = TRUE))
        query.words = query.words[query.words != ""]
        query.occ = matrix(rep(0, nrow(tfidf.table)), ncol = nrow(tfidf.table), nrow = 1)
        colnames(query.occ) = rownames(tfidf.table)
        for (i in query.words) {
            word.reg = paste0("^", i, "s?$")
            query.occ[1,str_detect(colnames(query.occ), regex(word.reg, ignore_case = TRUE))] = 1
        }
        
        # Compute cell type's score based on words match and words tfidf 
        ct.score = query.occ %*% tfidf.table
        db.match = colnames(ct.score)[ct.score[1,] >= threshold]
        
        # Try to find matches for a concatened form of the full query
        query.occ = matrix(rep(0, nrow(tfidf.table)), ncol = nrow(tfidf.table), nrow = 1)
        colnames(query.occ) = rownames(tfidf.table)
        query.reg = paste0(c("^", query.words, "s?$"), collapse = "")
        query.occ[1,str_detect(colnames(query.occ), regex(query.reg, ignore_case = TRUE))] = 1
        ct.score = query.occ %*% tfidf.table
        db.match.word = colnames(ct.score)[ct.score[1,] > 0.5]
        db.match = union(db.match, db.match.word)
        if (length(db.match) == 0) {
            print(paste0("No match for query: ", query))
        } 
    } 
    return(db.match)
}
