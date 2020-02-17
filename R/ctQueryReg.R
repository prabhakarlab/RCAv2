#' Function to parse from a character vector all character matching all keywords from a query.
#' 
#' @param query character or character vector.
#' @param db a character vector to parse query from. If remakedb=FALSE it should be a matrix with each extracted word of db as second column and with corresonding full character of db as first column.
#' @param remakedb a boolean value, default is TRUE. Put to FALSE (you should not do this) if you can provide a matrix as db as describe in db parameter.  
#' @return a character vector with all matches for one query.
#' @details 
#' This function will look for matches between a query and a character database using regular expression for text mining.
#' This function is not perfect and require some good practices:
#' - Take look at your database to design your query.
#' - Include both litteral and acronym form (if it exist) of a term in distinct query character.
#' - Avoid using plural forms in query.
#' - Correclty spell your query.
#' - Database will not always be clean, and concatenated words (for instance: forinstance)
#' will not be detected but with a specific corresponding concatened query.
#' @export
#' 
ctQueryReg = function(query, db ,remakedb=TRUE) {
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
    } else {
        new_db = db
    }
    # If user provides several queries for one db, apply this function for each query
    if (length(query) > 1) {
        db.match.list = apply(as.array(query), MARGIN = 1, FUN = "ct_query_reg", db=db, remakedb=FALSE)
        db.match = unlist(db.match.list)
        if (length(which(duplicated(db.match))) > 0) {
            db.match = db.match[-which(duplicated(db.match))]
        }
    } else {
        # Parse words in query
        query.words = gsub('([[:upper:]][[:lower:]])', ' \\1', query)
        query.words = as.vector(str_split(query.words, regex("[[:punct:] ]+"), simplify = TRUE))
        query.words = query.words[query.words != ""]
        db.match = c()
        # For each words in query, detect matching characters in db and 
        # Keep intersect with mactches of other words
        for (i in 1:length(query.words)) {
            query.reg = paste0("^", query.words[i], "s?$")
            db.match.word = new_db[str_detect(new_db[,2], regex(query.reg, ignore_case = TRUE)), 1]
            if (i == 1) {
                db.match = db.match.word
            } else {
                db.match = intersect(db.match, db.match.word)
            }
        }
        # Try to find matches for a concatened form of the full query
        query.reg = paste0(c("^", query.words, "S?$"), collapse = "")
        db.match.word = new_db[str_detect(new_db[,2], regex(query.reg, ignore_case = TRUE)), 1]
        db.match = union(db.match, db.match.word)
    } 
    return(db.match)
}
