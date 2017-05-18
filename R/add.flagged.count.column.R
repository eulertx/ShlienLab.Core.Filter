# this can be used within an apply statement to return
# the count of "FLAGGED" elements in a vector
add.flagged.count.column <- function(x) {
  count = 0
  for(i in 1:length(x)) {
    if (x[i] == 'FLAGGED') {
      count = count + 1
      }
    else {
      next()
      }
    }
  return(count)
  }
