#flip is a little function which is equal to flipud in matlab
# use recursive function to solve
flip = function(A) {
    n = dim(A)[1];
    if(is.null(n) || n==1) {
        return(A);
    }else {
        if(n>2) {
            return(rbind(A[n,],flip(A[2:(n-1),]),A[1,]));
        } else {
            return(rbind(A[n,],A[1,]));
        }
    }
}


